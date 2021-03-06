/*** rbm.c -- restricted boltzmann goodness
 *
 * Copyright (C) 2008-2013 Sebastian Freundt
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the author nor the names of any contributors
 *    may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ***/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <tgmath.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <setjmp.h>
#include <signal.h>
#include "maths.h"
#include "rand.h"
#include "nifty.h"

/* blas */
#if defined USE_BLAS
# include <mkl_cblas.h>
#endif	/* USE_BLAS */

#define DEFER_UPDATES

/* pick an implementation */
#if !defined SALAKHUTDINOV && !defined GEHLER
#define SALAKHUTDINOV	1
#endif	/* !SALAKHUTDINOV && !GEHLER */

#if defined __INTEL_COMPILER
# pragma warning (disable:1911)
# define auto		static
#endif	/* __INTEL_COMPILER */

#define PI		3.141592654f
#define E		2.718281828f

#if !defined NDEBUG
# define DEBUG(args...)	args
# define ni
#else  /* NDEBUG */
# define DEBUG(args...)
# define ni
#endif	/* !NDEBUG */


/* my alibi blas */
#if !defined USE_BLAS
typedef long int MKL_INT;
#endif	/* !USE_BLAS */

#if defined __SSE__
static ni float
sse_sdot(const MKL_INT N, const float *X, const float *Y)
{
#define f4z		sizeof(float) * 4U
	typedef float f4 __attribute__((aligned(f4z), vector_size(f4z)));
	const f4 *X4 = (const f4*)X;
	const f4 *Y4 = (const f4*)Y;
	union {
		f4 v;
		float f[4U];
	} sum = {};

	for (MKL_INT i = 0; i < N / 4U; i++, X++, Y++) {
		sum.v += *X4 * *Y4;
	}
	return sum.f[0U] + sum.f[1U] + sum.f[2U] + sum.f[3U];
}
#endif	/* __SSE__ */

static ni float
drb_sdot11(const MKL_INT N, const float *X, const float *Y)
{
	float sum = 0.f;

#if defined __SSE__
	if (!(N % 4U)) {
		/* ah, we can use the sse version */
		return sse_sdot(N, X, Y);
	}
#endif	/* __SSE__ */
	/* otherwise proceed with the bog standard procedure */
	for (MKL_INT i = 0; i < N; i++, X++, Y++) {
		sum += *X * *Y;
	}
	return sum;
}

static ni float*
tr(const float *w, const MKL_INT m, const MKL_INT n)
{
	float *res = malloc(m * n * sizeof(*res));

#define wtr(i, j)	res[i * m + j]
#define w(i, j)		w[i * n + j]
	for (MKL_INT i = 0; i < m; i++) {
		for (MKL_INT j = 0; j < n; j++) {
			wtr(j, i) = w(i, j);
		}
	}
#undef w
#undef wtr
	return res;
}


/* mmapping, adapted from fops.h */
typedef struct glodf_s glodf_t;
typedef struct glodfn_s glodfn_t;

struct glodf_s {
	size_t z;
	void *d;
};

struct glodfn_s {
	int fd;
	struct glodf_s fb;
};

static inline glodf_t
mmap_fd(int fd, size_t fz, int prot, int flags)
{
	void *p;

	p = mmap(NULL, fz, prot, flags, fd, 0);
	if (p == MAP_FAILED) {
		return (glodf_t){.z = 0U, .d = NULL};
	}
	return (glodf_t){.z = fz, .d = p};
}

static inline int
munmap_fd(glodf_t map)
{
	return munmap(map.d, map.z);
}

static glodfn_t
mmap_fn(const char *fn, int flags)
{
	const int fl = MAP_SHARED;
	const int pr = PROT_READ | ((flags & O_RDWR) ? PROT_WRITE : 0);
	struct stat st;
	glodfn_t res;

	if ((res.fd = open(fn, flags)) < 0) {
		;
	} else if (fstat(res.fd, &st) < 0) {
		res.fb = (glodf_t){.z = 0U, .d = NULL};
		goto clo;
	} else if ((res.fb = mmap_fd(res.fd, st.st_size, pr, fl)).d == NULL) {
	clo:
		close(res.fd);
		res.fd = -1;
	}
	return res;
}

static int
munmap_fn(glodfn_t f)
{
	int rc = 0;

	if (f.fb.d != NULL) {
		rc += munmap_fd(f.fb);
	}
	if (f.fd >= 0) {
		rc += close(f.fd);
	}
	return rc;
}

static glodfn_t
mremap_fn(glodfn_t f, int nu_flags)
{
	const int nu_prot = PROT_READ;

	if (UNLIKELY(f.fb.d == NULL)) {
		return (glodfn_t){.fd = -1};
	}
	/* first munmap the old buffer */
	munmap_fd(f.fb);
	f.fb = mmap_fd(f.fd, f.fb.z, nu_prot, nu_flags);
	return f;
}


/* custom machine */
typedef struct dl_rbm_s *dl_rbm_t;

struct dl_rbm_s {
	size_t nvis;
	size_t nhid;
	float *vbias;
	float *hbias;
	float *w;

	void *priv;
};

struct dl_spec_s {
	size_t nvis;
	size_t nhid;
};

struct dl_file_s {
	uint8_t magic[4U];
	uint8_t flags[4U];
	struct dl_spec_s sp;

	/* offset to the beginning */
	size_t off;
	float data[];
};

struct dl_rbm_priv_s {
	glodfn_t f;
	/* transpose of w */
	float *wtr;
};

static dl_rbm_t
pump(const char *file, int flags)
{
	static struct dl_rbm_s res;
	static struct dl_rbm_priv_s p[1];

	if (UNLIKELY((p->f = mmap_fn(file, flags)).fd < 0)) {
		goto out;
	} else if (UNLIKELY(p->f.fb.z < sizeof(struct dl_file_s))) {
		goto out;
	}
	with (struct dl_file_s *fl = p->f.fb.d) {
		float_t *dp = fl->data + fl->off;

		res.nvis = fl->sp.nvis;
		res.vbias = dp;
		dp += fl->sp.nvis;

		res.nhid = fl->sp.nhid;
		res.hbias = dp;
		dp += fl->sp.nhid;

		res.w = dp;

		p->wtr = tr(res.w, fl->sp.nvis, fl->sp.nhid);
	}
	res.priv = p;
	return &res;
out:
	/* and out are we */
	(void)munmap_fn(p->f);
	return NULL;
}

static int
dump(dl_rbm_t m)
{
	struct dl_rbm_priv_s *p = m->priv;

	if (UNLIKELY(m == NULL)) {
		return 0;
	}
	return munmap_fn(p->f);
}

static int
resz(dl_rbm_t m, struct dl_spec_s nu)
{
/* shrink or expand the machine in M according to dimensions in NU */
	const int pr = PROT_READ | PROT_WRITE;
	const int fl = MAP_SHARED;
	struct dl_rbm_priv_s *p = m->priv;
	glodfn_t ol_f;
	glodfn_t nu_f;
	size_t z;
	size_t fz;
	int res = 0;

	/* remap the current guy */
	ol_f = mremap_fn(p->f, MAP_PRIVATE);

	/* compute new file size */
	z = nu.nvis + nu.nhid + nu.nvis * nu.nhid;
	ftruncate(ol_f.fd, fz = z * sizeof(float) + sizeof(struct dl_file_s));

	if (UNLIKELY((nu_f.fb = mmap_fd(ol_f.fd, fz, pr, fl)).d == NULL)) {
		res = -1;
		goto out;
	}
	nu_f.fd = ol_f.fd;

	/* just copy the file header over */
	with (struct dl_file_s *fp = nu_f.fb.d, *op = ol_f.fb.d) {
		const float vnois = .1f;
		const float hnois = .01f;
		const float wnois = 1.f / (nu.nvis * nu.nhid);
		/* sp is our source pointer (in the private map) */
		const float_t *sp = op->data + op->off;
		/* dp is our target pointer in the truncated file */
		float_t *restrict dp = fp->data + fp->off;

		/* vbiasses */
		memcpy(m->vbias = dp, sp, op->sp.nvis * sizeof(*sp));
		dp += m->nvis = nu.nvis;
		sp += op->sp.nvis;
		/* wobble */
		for (size_t i = op->sp.nvis; i < m->nvis; i++) {
			const float x = dr_rand_uni();
			m->vbias[i] = log(vnois * x);
		}

		/* hbiasses */
		memcpy(m->hbias = dp, sp, op->sp.nhid * sizeof(*sp));
		dp += m->nhid = nu.nhid;
		sp += op->sp.nhid;
		/* wobble */
		for (size_t j = op->sp.nhid; j < m->nhid; j++) {
			const float x = dr_rand_norm();
			m->hbias[j] = hnois * x;
		}

		/* weight matrix */
		memcpy(m->w = dp, sp, op->sp.nvis * op->sp.nhid * sizeof(*sp));
		/* wobble */
		for (size_t k = op->sp.nvis * op->sp.nhid;
		     k < m->nvis * m->nhid; k++) {
			const float x = dr_rand_norm();
			m->w[k] = wnois * x;
		}

		/* recalc the transposed of w */
		free(p->wtr);
		p->wtr = tr(m->w, nu.nvis, nu.nhid);

		/* assign the new size */
		fp->sp = nu;
	}
	p->f = nu_f;

out:
	munmap_fd(ol_f.fb);
	return res;
}

static dl_rbm_t
crea(const char *file, struct dl_spec_s sp)
{
	static glodfn_t f;
	const int pr = PROT_READ | PROT_WRITE;
	const int fl = MAP_SHARED;
	size_t fz;
	int fd;
	dl_rbm_t res;

	if (UNLIKELY((fd = open(file, O_CREAT | O_RDWR | O_TRUNC, 0666)) < 0)) {
		return NULL;
	}
	/* compute total file size */
	ftruncate(fd, fz = sizeof(struct dl_file_s));

	if (UNLIKELY((f.fb = mmap_fd(f.fd = fd, fz, pr, fl)).d == NULL)) {
		goto out;
	}
	/* normally we'd fill in titbits like our magic string and flags */
	;
	munmap_fn(f);

	/* now let pump() and resz() do the rest */
	if ((res = pump(file, O_RDWR)) == NULL) {
		/* nawww */
		goto out;
	} else if (resz(res, sp) < 0) {
		/* shame */
		goto out;
	}
	/* all's good */
	return res;

out:
	close(fd);
	return NULL;
}


/* sparse integer vectors */
typedef struct spsc_s spsc_t;
typedef struct spsv_s spsv_t;

struct spsc_s {
	size_t i;
	uint8_t v;
};

struct spsv_s {
	size_t z;
	spsc_t *v;
};

static spsv_t
read_tf(const int fd)
{
	static char *line;
	static size_t llen;
	static struct spsc_s *spsv;
	static size_t spsz;
	ssize_t nrd;
	size_t i = 0U;

	if (UNLIKELY(fd < 0)) {
		if (line != NULL) {
			free(line);
		}
		if (spsv != NULL) {
			free(spsv);
		}
		return (spsv_t){.z = 0U, .v = NULL};
	}

	/* now then */
	while ((nrd = getline(&line, &llen, stdin)) > 0) {
		char *p;
		long unsigned int v;
		long unsigned int c;

		/* check for form feeds, and maybe yield */
		if (*line == '\f') {
			break;
		}
		/* read the term id */
		v = strtoul(line, &p, 0);
		if (*p++ != '\t') {
			continue;
		}
		/* read the count */
		c = strtoul(p, &p, 0);
		if (*p++ != '\n') {
			continue;
		}

		if (UNLIKELY(i >= spsz)) {
			/* extend vector */
			size_t nu = spsz + 256U;
			spsv = realloc(spsv, nu * sizeof(*spsv));
			spsz = nu;
		}
		/* assign index/value pair */
		spsv[i].i = v;
		spsv[i].v = c;
		i++;
	}
	return (spsv_t){.z = i, .v = spsv};
}

static size_t
popul_sv(float *restrict x, size_t z, const spsv_t sv)
{
/* Populate the bottom visible layer X (hopefully large enough)
 * with values from sparse vector SV.
 * Return the total number of words. */
	size_t res = 0U;

	memset(x, 0, z * sizeof(*x));
	for (size_t j = 0; j < sv.z; j++) {
		size_t i = sv.v[j].i;
		uint8_t c = sv.v[j].v;

		if (UNLIKELY(i >= z)) {
			fprintf(stderr, "\
not populating entry %zu, machine's network too small\n", i);
			continue;
		}

		res += c;
		x[i] = (float)(int)c;
	}
	return res;
}


#if !defined NDEBUG
static void
dump_layer(const char *pre, const float *x, size_t z)
{
	float minx = INFINITY;
	float maxx = -INFINITY;

	for (size_t i = 0; i < z; i++) {
		if (x[i] < minx) {
			minx = x[i];
		}
		if (x[i] > maxx) {
			maxx = x[i];
		}
	}
	printf("%s (%.6g  %.6g)\n", pre, minx, maxx);
	return;
}

static size_t
count_layer(const float *x, size_t z)
{
	size_t sum = 0U;

	for (size_t i = 0; i < z; i++) {
		if (x[i] > 0.f) {
			sum++;
		}
	}
	return sum;
}

static float
integ_layer(const float *x, size_t z)
{
	float sum = 0U;

	for (size_t i = 0; i < z; i++) {
		if (x[i] > 0.f) {
			sum += x[i];
		}
	}
	return sum;
}
#endif	/* !NDEBUG */


/* propagation, gibbs sampling and learning */
/* global parameters */
#if defined SALAKHUTDINOV
static size_t N;
#endif	/* SALAKHUTDINOV */

static ni int
prop_up(float *restrict h, dl_rbm_t m, const float vis[static m->nvis])
{
/* propagate visible units activation upwards to the hidden units (recon) */
	const size_t nvis = m->nvis;
	const size_t nhid = m->nhid;
	const struct dl_rbm_priv_s *p = m->priv;
	const float *wtr = p->wtr;
	const float *b = m->hbias;

#define w(j)		(wtr + j * nvis)
	for (size_t j = 0; j < nhid; j++) {
		h[j] = b[j] + drb_sdot11(nvis, w(j), vis);
	}
#undef w
	return 0;
}

static ni int
expt_hid(float *restrict h, dl_rbm_t m, const float hid[static m->nhid])
{
	const size_t nhid = m->nhid;

	DEBUG(dump_layer("Ha", hid, nhid));

	for (size_t j = 0; j < nhid; j++) {
		h[j] = sigma(hid[j]);
	}
	return 0;
}

static ni int
smpl_hid(float *restrict h, dl_rbm_t m, const float hid[static m->nhid])
{
/* infer hidden unit states given vis(ible units) */
	const size_t nhid = m->nhid;

	DEBUG(dump_layer("He", hid, nhid));

	for (size_t j = 0; j < nhid; j++) {
		/* just flip a coin */
		h[j] = dr_rand_binom1(hid[j]);
	}

	DEBUG(dump_layer("Hs", h, nhid));
	return 0;
}

static ni int
prop_down(float *restrict v, dl_rbm_t m, const float hid[static m->nhid])
{
/* propagate hidden units activation downwards to the visible units */
	const size_t nvis = m->nvis;
	const size_t nhid = m->nhid;
	const float *w = m->w;
	const float *b = m->vbias;

#define w(i)		(w + i * nhid)
	for (size_t i = 0; i < nvis; i++) {
		v[i] = b[i] + drb_sdot11(nhid, w(i), hid);
	}
#undef w
	return 0;
}

static ni int
expt_vis(float *restrict v, dl_rbm_t m, const float vis[static m->nvis])
{
	const size_t nvis = m->nvis;

	DEBUG(dump_layer("Va", vis, nvis));

#if defined SALAKHUTDINOV
	softmax(v, vis, nvis);
	with (const float scal = (float)N) {
		for (size_t i = 0; i < nvis; i++) {
			v[i] *= scal;
		}
	}
#elif defined GEHLER
	for (size_t i = 0; i < nvis; i++) {
		v[i] = exp(vis[i]);
	}
#elif defined BINOM_INPUT
	for (size_t i = 0; i < nvis; i++) {
		v[i] = sigma(vis[i]);
	}
#endif	/* impls */
	return 0;
}

static ni int
smpl_vis(float *restrict v, dl_rbm_t m, const float vis[static m->nvis])
{
/* infer visible unit states given hid(den units) */
	const size_t nvis = m->nvis;

	DEBUG(dump_layer("Ve", vis, nvis));

	/* vis is expected to contain the lambda values */
	for (size_t i = 0; i < nvis; i++) {
#if !defined BINOM_INPUT
		v[i] = dr_rand_poiss(vis[i]);
#else  /* BINOM_INPUT */
		v[i] = dr_rand_binom1(vis[i]);
#endif	/* !BINOM_INPUT */
	}

	DEBUG(dump_layer("Vs", vis, nvis));
	return 0;
}


/* training and classifying modes */
typedef struct drbctx_s *drbctx_t;

struct drbctx_s {
	dl_rbm_t m;

	/* some scratch vectors for the gibbs sampling */
	float *vo;
	float *ho;
	float *vr;
	float *hr;

	/* difference vectors */
	float *dh;
	float *dv;
	float *dw;
};

static const float eta = 0.02f;
static const float mom = 0.9f;
static const float dec = 0.f;

static void
init_drbctx(struct drbctx_s *restrict tgt, dl_rbm_t m)
{
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;

	tgt->m = m;

	/* initialise the scratch vectors */
	tgt->vo = calloc(nv, sizeof(*tgt->vo));
	tgt->vr = calloc(nv, sizeof(*tgt->vr));
	tgt->ho = calloc(nh, sizeof(*tgt->ho));
	tgt->hr = calloc(nh, sizeof(*tgt->hr));

	tgt->dw = calloc(nh * nv, sizeof(*tgt->dw));
	tgt->dh = calloc(nh, sizeof(*tgt->dh));
	tgt->dv = calloc(nv, sizeof(*tgt->dv));
	return;
}

static void
fini_drbctx(struct drbctx_s *tgt)
{
	tgt->m = NULL;

	free(tgt->vo);
	free(tgt->vr);
	free(tgt->ho);
	free(tgt->hr);

	free(tgt->dw);
	free(tgt->dv);
	free(tgt->dh);
	return;
}

static void
rset_drbctx(struct drbctx_s *tgt)
{
	const size_t nv = tgt->m->nvis;
	const size_t nh = tgt->m->nhid;

	memset(tgt->dw, 0, nv * nh);
	memset(tgt->dv, 0, nv);
	memset(tgt->dh, 0, nh);
	return;
}

#if defined __SSE__ && defined DEFER_UPDATES && 0
/* this version is slower than the sequential one */
static ni void
sse_update_w(drbctx_t ctx)
{
#define f4z		sizeof(float) * 4U
	typedef float f4 __attribute__((aligned(f4z), vector_size(f4z)));
	dl_rbm_t m = ctx->m;
	const float *vo = ctx->vo;
	const float *ho = ctx->ho;
	const float *vr = ctx->vr;
	const float *hr = ctx->hr;
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;
	float *restrict dw = ctx->dw;
	const float *w = m->w;
#if !defined NDEBUG
	float mind = INFINITY;
	float maxd = -INFINITY;
#endif	/* !NDEBUG */

#define w(i, j)		(w + i * nh + j)
#define dw(i, j)	(dw + i * nh + j)

	/* bang <v_i h_j> into weights */
	for (size_t i = 0; i < nv; i++) {
		for (size_t j = 0; j < nh; j += 4U) {
			const f4 vo4 = {vo[i], vo[i], vo[i], vo[i]};
			const f4 ho4 = *(const f4*)(ho + j);
			const f4 vr4 = {vr[i], vr[i], vr[i], vr[i]};
			const f4 hr4 = *(const f4*)(hr + j);
			const f4 vho4 = vo4 * ho4;
			const f4 vhr4 = vr4 * hr4;
			f4 d = vho4 - vhr4;
			static const f4 dec4 = {0.f, 0.f, 0.f, 0.f};
			static const f4 eta4 = {0.02f, 0.02f, 0.02f, 0.02f};
			static const f4 mom4 = {0.9f, 0.9f, 0.9f, 0.9f};

			/* decay */
			d -= dec4 * *(const f4*)w(i, j);
			/* learning rate */
			d *= eta4;
			/* momentum term */
			d += mom4 * *(const f4*)dw(i, j);

#if !defined NDEBUG && 0
			if (d < mind) {
				mind = d;
			}
			if (d > maxd) {
				maxd = d;
			}
#endif	/* !NDEBUG */
			*(f4*)dw(i, j) = d;
		}
	}
#if !defined NDEBUG
	printf("dw (%.6g  %.6g)\n", mind, maxd);
#endif	/* !NDEBUG */
#undef w
#undef dw
	return;
}
#endif	/* __SSE__ */

static ni void
update_w(drbctx_t ctx)
{
	dl_rbm_t m = ctx->m;
	const struct dl_rbm_priv_s *p = m->priv;
	const float *vo = ctx->vo;
	const float *ho = ctx->ho;
	const float *vr = ctx->vr;
	const float *hr = ctx->hr;
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;
	float *restrict dw = ctx->dw;
#if defined DEFER_UPDATES
	const float *w = m->w;
	const float *UNUSED(wtr) = p->wtr;
#else  /* !DEFER_UPDATES */
	float *restrict w = m->w;
	float *restrict wtr = p->wtr;
#endif	/* DEFER_UPDATES */
#if !defined NDEBUG
	float mind = INFINITY;
	float maxd = -INFINITY;
#endif	/* !NDEBUG */

#define w(i, j)		w[i * nh + j]
#define wtr(i, j)	wtr[i * nv + j]
#define dw(i, j)	dw[i * nh + j]

	/* bang <v_i h_j> into weights */
	for (size_t i = 0; i < nv; i++) {
		for (size_t j = 0; j < nh; j++) {
			float vho = vo[i] * ho[j];
			float vhr = vr[i] * hr[j];
			float d = vho - vhr;

			/* decay */
			d -= dec * w(i, j);
			/* learning rate */
			d *= eta;
			/* momentum term */
			d += mom * dw(i, j);

#if !defined NDEBUG
			if (d < mind) {
				mind = d;
			}
			if (d > maxd) {
				maxd = d;
			}
#endif	/* !NDEBUG */
#if !defined DEFER_UPDATES
			wtr(j, i) += d;
			w(i, j) +=
#endif	/* !DEFER_UPDATES */
				dw(i, j) = d;
		}
	}
#if !defined NDEBUG
	printf("dw (%.6g  %.6g)\n", mind, maxd);
#endif	/* !NDEBUG */
#undef w
#undef wtr
#undef dw
	return;
}

static ni void
update_b(drbctx_t ctx)
{
	dl_rbm_t m = ctx->m;
	const float *vo = ctx->vo;
	const float *ho = ctx->ho;
	const float *vr = ctx->vr;
	const float *hr = ctx->hr;
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;

	/* the real update routine is here, called twice, for vb and hb */
	auto void __upd(
		float *restrict b, float *restrict db,
		const float *o, const float *r, const size_t z)
	{
#if !defined NDEBUG
		float mind = INFINITY;
		float maxd = -INFINITY;
#endif	/* !NDEBUG */

		for (size_t i = 0; i < z; i++) {
			float d = o[i] - r[i];

			/* decay */
			d -= dec * b[i];
			/* learning rate */
			d *= eta;
			/* momentum */
			d += mom * db[i];

#if !defined NDEBUG
			if (d < mind) {
				mind = d;
			}
			if (d > maxd) {
				maxd = d;
			}
#endif	/* !NDEBUG */
#if !defined DEFER_UPDATES
			b[i] +=
#endif	/* !DEFER_UPDATES */
				db[i] = d;
		}
#if !defined NDEBUG
		printf("db (%.6g  %.6g)\n", mind, maxd);
#endif	/* !NDEBUG */
	}

	/* bias update */
	__upd(m->vbias, ctx->dv, vo, vr, nv);
	__upd(m->hbias, ctx->dh, ho, hr, nh);
	return;
}

static ni void
final_update_w(drbctx_t ctx)
{
/* finalise the weight update */
#if defined DEFER_UPDATES
	const size_t nv = ctx->m->nvis;
	const size_t nh = ctx->m->nhid;

#define w(i, j)		ctx->m->w[i * nh + j]
#define dw(i, j)	ctx->dw[i * nh + j]

	/* now really bang <v_i h_j> into weights */
	for (size_t i = 0; i < nv; i++) {
		for (size_t j = 0; j < nh; j++) {
			w(i, j) += dw(i, j);
		}
	}
#undef w
#undef dw
#else  /* !DEFER_UPDATES */
	ctx = ctx;
#endif	/* DEFER_UPDATES */
	return;
}

static ni void
final_update_b(drbctx_t ctx)
{
/* finalise the weight update */
#if defined DEFER_UPDATES
	const size_t nv = ctx->m->nvis;
	const size_t nh = ctx->m->nhid;

#define v(i)		ctx->m->vbias[i]
#define dv(i)		ctx->dv[i]

	/* now really bang bias updates into the biasses */
	for (size_t i = 0; i < nv; i++) {
		v(i) += dv(i);
	}
#undef v
#undef dv

#define h(i)		ctx->m->hbias[i]
#define dh(i)		ctx->dh[i]
	for (size_t j = 0; j < nh; j++) {
		h(j) += dh(j);
	}
#else  /* !DEFER_UPDATES */
	ctx = ctx;
#endif	/* DEFER_UPDATES */
	return;
}

static void
train(drbctx_t ctx, struct spsv_s sv)
{
	const dl_rbm_t m = ctx->m;
	float *restrict vo = ctx->vo;
	float *restrict ho = ctx->ho;
	float *restrict vr = ctx->vr;
	float *restrict hr = ctx->hr;
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;

	DEBUG(float *hs = calloc(nh, sizeof(*hs)));

	/* populate from input */
#if defined SALAKHUTDINOV
	N = popul_sv(vo, nv, sv);
#else  /* !SALAKHUTDINOV */
	(void)popul_sv(vo, nv, sv);
#endif	/* SALAKHUTDINOV */

	/* vh gibbs */
	prop_up(ho, m, vo);
	expt_hid(ho, m, ho);
	/* don't sample into ho, use hr instead, we want the activations */
	smpl_hid(hr, m, ho);
	DEBUG(size_t nho = count_layer(hr, nh));

	/* hv gibbs */
	prop_down(vr, m, hr);
	expt_vis(vr, m, vr);
	smpl_vis(vr, m, vr);
	/* vh gibbs */
	prop_up(hr, m, vr);
	expt_hid(hr, m, hr);

	DEBUG(
		smpl_hid(hs, m, hr);
		size_t nhr = count_layer(hs, nh);
		);

#if !defined NDEBUG
	size_t nso = count_layer(vo, nv);
	size_t nsr = count_layer(vr, nv);
	float Nso = integ_layer(vo, nv);
	float Nsr = integ_layer(vr, nv);
	printf("|vo| %zu  |vr| %zu  Nvo %.6g  Nvr %.6g\n", nso, nsr, Nso, Nsr);
	printf("|ho| %zu  |hr| %zu\n", nho, nhr);
#endif	/* !NDEBUG */

	/* we won't sample the h reconstruction as we want to use the
	 * the activations directly */
	update_w(ctx);
	update_b(ctx);

#if !defined NDEBUG
	dump_layer("h", m->hbias, nh);
	dump_layer("v", m->vbias, nv);
	dump_layer("w", m->w, nh * nv);
#endif	/* !NDEBUG */

	DEBUG(free(hs));
	return;
}

static void
prop(drbctx_t ctx, spsv_t sv, int smplp)
{
#define m	ctx->m
#define vo	ctx->vo
#define ho	ctx->ho
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;

	/* populate from input */
#if defined SALAKHUTDINOV
	N = popul_sv(vo, nv, sv);
#else  /* !SALAKHUTDINOV */
	(void)popul_sv(vo, nv, sv);
#endif	/* SALAKHUTDINOV */

	/* vh gibbs */
	prop_up(ho, m, vo);
	expt_hid(ho, m, ho);
	if (!smplp) {
		for (size_t i = 0; i < nh; i++) {
			printf("%g\n", ho[i]);
		}
	} else {
		/* oh, madame wants sampling as well */
		size_t nsmpl = 0U;

		smpl_hid(ho, m, ho);

		for (size_t i = 0U; i < nh; i++) {
			uint8_t hi = (uint8_t)(int)ho[i];

			if (UNLIKELY(hi)) {
				printf("%zu\t%u\n", i, (unsigned int)hi);
				nsmpl++;
			}
		}
		if (LIKELY(nsmpl)) {
			puts("\f");
		}
	}
	return;
#undef m
#undef vo
#undef ho
}

static int
check(dl_rbm_t m)
{
	int res = 0;

	with (size_t nv = m->nvis, nh = m->nhid) {
		for (size_t i = 0; i < nv; i++) {
			if (UNLIKELY(isnan(m->vbias[i]))) {
				printf("VBIAS[%zu] <- NAN\n", i);
				res = 1;
			}
		}
		for (size_t j = 0; j < nh; j++) {
			if (UNLIKELY(isnan(m->hbias[j]))) {
				printf("HBIAS[%zu] <- NAN\n", j);
				res = 1;
			}
		}
		for (size_t i = 0; i < nv; i++) {
			for (size_t j = 0; j < nh; j++) {
				if (UNLIKELY(isnan(m->w[i * nh + j]))) {
					printf("W[%zu,%zu] <- NAN\n", i, j);
					res = 1;
				}
			}
		}
	}
	return res;
}


#if defined __INTEL_COMPILER
# pragma warning (disable:593)
# pragma warning (disable:181)
#endif	/* __INTEL_COMPILER */
#include "rbm.xh"
#include "rbm.x"
#if defined __INTEL_COMPILER
# pragma warning (default:593)
# pragma warning (default:181)
#endif	/* __INTEL_COMPILER */

static int
cmd_init(struct glod_args_info argi[static 1])
{
/* shell return codes, 0 success, 1 failure */
	const char *file = argi->inputs[1U];
	struct dl_spec_s dim;
	dl_rbm_t m;
	int res = 0;

	/* parse dimens */
	with (char *chk = argi->dimen_arg) {
		dim.nvis = strtoul(chk, &chk, 0);
		if ((*chk++ | 0x20) != 'x') {
			res = 1;
			goto out;
		}
		dim.nhid = strtoul(chk, &chk, 0);
		if (*chk) {
			res = 1;
			goto out;
		}
	}

	/* just create (or resize) the machine */
	if (argi->inputs_num < 2) {
		fputs("no machine file given\n", stderr);
		res = 1;
	} else if (!argi->resize_given && (m = crea(file, dim)) == NULL) {
		fprintf(stderr, "error creating machine file `%s'\n", file);
		res = 1;
	} else if (argi->resize_given && (m = pump(file, O_RDWR)) == NULL) {
		fprintf(stderr, "error loading machine file `%s'\n", file);
		res = 1;
	} else if (argi->resize_given && resz(m, dim) < 0) {
		/* actually do resize now (in --resize mode) */
		res = 1;
	} else {
		/* for both cases, write the file to disk */
		res = dump(m);
	}
out:
	return res;
}

static int
cmd_train(struct glod_args_info argi[static 1])
{
	static jmp_buf jb;
	const char *file = argi->inputs[1U];
	dl_rbm_t m = NULL;
	int res = 0;

	if (argi->inputs_num < 2) {
		fputs("no machine file given\n", stderr);
		res = 1;

	} else if (UNLIKELY((m = pump(file, O_RDWR)) == NULL)) {
		/* reading the machine file failed */
		fprintf(stderr, "error opening machine file `%s'\n", file);
		res = 1;

	} else if (setjmp(jb)) {
		/* C-c handler */
		goto train_xit;

	} else {
		/* all clear */
		static struct drbctx_s ctx[1];
		const int fd = STDIN_FILENO;
		const size_t batchz = argi->batch_size_arg;
		size_t i = 0;

		/* set up C-c handling (dirtee) */
		auto __attribute__((noreturn)) void si_train(int UNUSED(sig))
		{
			longjmp(jb, 1U);
		}
		signal(SIGINT, si_train);

		init_rand();
		init_drbctx(ctx, m);

		for (spsv_t sv; (sv = read_tf(fd)).z; train(ctx, sv)) {
			if (++i == batchz) {
				/* update weights and biasses */
				final_update_w(ctx);
				final_update_b(ctx);
				rset_drbctx(ctx);
				i = 0U;
			}
		}
	train_xit:
		/* also hopped to by the signal handler */
		final_update_w(ctx);
		final_update_b(ctx);

		/* just to deinitialise resources */
		(void)read_tf(-1);
		fini_drbctx(ctx);
		dump(m);
		deinit_rand();
	}
	return res;
}

static int
cmd_prop(struct glod_args_info argi[static 1])
{
	static jmp_buf jb;
	const char *file = argi->inputs[1U];
	dl_rbm_t m = NULL;
	int res = 0;

	if (argi->inputs_num < 2) {
		fputs("no machine file given\n", stderr);
		res = 1;

	} else if (UNLIKELY((m = pump(file, O_RDONLY)) == NULL)) {
		/* reading the machine file failed */
		fprintf(stderr, "error opening machine file `%s'\n", file);
		res = 1;

	} else if (setjmp(jb)) {
		/* C-c handler */
		goto prop_xit;

	} else {
		/* all clear */
		static struct drbctx_s ctx[1];
		const int fd = STDIN_FILENO;
		const int smplp = argi->sample_given;

		/* set up the C-c handler */
		auto __attribute__((noreturn)) void si_prop(int UNUSED(sig))
		{
			longjmp(jb, 2U);
		}
		signal(SIGINT, si_prop);

		init_rand();
		init_drbctx(ctx, m);

		for (spsv_t sv; (sv = read_tf(fd)).z; prop(ctx, sv, smplp));

	prop_xit:
		/* just to deinitialise resources */
		(void)read_tf(-1);
		fini_drbctx(ctx);
		dump(m);
		deinit_rand();
	}

	return res;
}

static int
cmd_info(struct glod_args_info argi[static 1])
{
	int res = 0;

	for (unsigned int i = 1; i < argi->inputs_num; i++) {
		const char *f = argi->inputs[i];
		dl_rbm_t m;

		if ((m = pump(f, O_RDONLY)) == NULL) {
			fprintf(stderr, "error opening machine file `%s'\n", f);
			res = 1;
		}

		/* just a general overview */
		printf("%s\t%zux%zu\tpoiss->binary\n", f, m->nvis, m->nhid);

		/* and close the resources assoc'd with m again */
		dump(m);
	}
	return res;
}


int
main(int argc, char *argv[])
{
	struct glod_args_info argi[1];
	int res = 0;

	if (glod_parser(argc, argv, argi)) {
		res = 1;
		goto out;
	} else if (argi->inputs_num < 1) {
		glod_parser_print_help();
		res = 1;
		goto out;
	}

	/* check the command */
	with (const char *cmd = argi->inputs[0]) {
		if (!strcmp(cmd, "train")) {
			res = cmd_train(argi);

		} else if (!strcmp(cmd, "prop")) {
			res = cmd_prop(argi);

		} else if (!strcmp(cmd, "init")) {
			res = cmd_init(argi);

		} else if (!strcmp(cmd, "info")) {
			res = cmd_info(argi);

		} else {
			/* otherwise print help and bugger off */
			glod_parser_print_help();
			res = 1;
		}
	}

out:
	glod_parser_free(argi);
	return res;
}

/* rbm.c ends here */
