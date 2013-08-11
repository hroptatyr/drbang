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
#include "rand.h"
#include "nifty.h"

/* blas */
#if defined USE_BLAS
# include <mkl_cblas.h>
#endif	/* USE_BLAS */

#define PREFER_NUMERICAL_STABILITY_OVER_SPEED
#define DEFER_UPDATES

/* pick an implementation */
#if !defined SALAKHUTDINOV && !defined GEHLER
#define SALAKHUTDINOV	1
#endif	/* !SALAKHUTDINOV && !GEHLER */

#if defined __INTEL_COMPILER
# pragma warning (disable:1911)
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


static float
factorialf(uint8_t n)
{
	static const float table[] = {
		1., 1., 2., 6., 24., 120., 720., 5040.,
		40320., 362880., 3628800., 39916800.,
		479001600., 6227020800., 87178291200., 1307674368000.,
	};
	float res;

	if (LIKELY(n < countof(table))) {
		return table[n];
	}

	/* otherwise proceed from 16! */
	res = table[15U];
	n -= 15U;
	for (float x = 16.; n > 0; n--, x += 1.0) {
		res *= x;
	}
	return res;
}

static double
factorial(uint8_t n)
{
	static const double table[] = {
		1., 1., 2., 6., 24., 120., 720., 5040.,
		40320., 362880., 3628800., 39916800.,
		479001600., 6227020800., 87178291200., 1307674368000.,
	};
	double res;

	if (LIKELY(n < countof(table))) {
		return table[n];
	}

	/* otherwise proceed from 16! */
	res = table[15U];
	n -= 15U;
	for (double x = 16.; n > 0; n--, x += 1.0) {
		res *= x;
	}
	return res;
}

static long double
factoriall(uint8_t n)
{
	static const long double table[] = {
		1.l, 1.l, 2.l, 6.l, 24.l, 120.l, 720.l, 5040.l,
		40320.l, 362880.l, 3628800.l, 39916800.l,
		479001600.l, 6227020800.l, 87178291200.l, 1307674368000.l,
	};
	long double res;

	if (LIKELY(n < countof(table))) {
		return table[n];
	}

	/* otherwise proceed from 16! */
	res = table[15U];
	n -= 15U;
	for (long double x = 16.; n > 0; n--, x += 1.0) {
		res *= x;
	}
	return res;
}

static __attribute__((unused)) float
poissf(float lambda, uint8_t n)
{
	float res;

	res = exp(-lambda);
	for (uint8_t i = n; i > 0; i--) {
		res *= lambda;
	}
	res /= factorialf(n);
	return res;
}

static __attribute__((unused)) double
poiss(double lambda, uint8_t n)
{
	double res;

	res = exp(-lambda);
	for (uint8_t i = n; i > 0; i--) {
		res *= lambda;
	}
	res /= factorial(n);
	return res;
}

static __attribute__((unused)) long double
poissl(long double lambda, uint8_t n)
{
	long double res;

	res = exp(-lambda);
	for (uint8_t i = n; i > 0; i--) {
		res *= lambda;
	}
	res /= factoriall(n);
	return res;
}

static float
sigmaf(float x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1.f + tanh(x / 2.f)) / 2.f;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1.f / (1.f + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

static double
sigma(double x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1. + tanh(x / 2.)) / 2.;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1. / (1. + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

static long double
sigmal(long double x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1.L + tanh(x / 2.L)) / 2.L;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1.L / (1.L + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

/* my own tgmaths */
#define poiss(x, n)	__TGMATH_BINARY_FIRST_REAL_ONLY(x, n, poiss)
#define sigma(x)	__TGMATH_UNARY_REAL_ONLY(x, sigma)


/* my alibi blas */
#if !defined USE_BLAS
typedef long int MKL_INT;

static ni float
cblas_sdot(
	const MKL_INT N,
	const float *X, const MKL_INT incX,
	const float *Y, const MKL_INT incY)
{
	float sum = 0.f;

	for (MKL_INT i = 0; i < N; i++, X += incX, Y += incY) {
		sum += *X * *Y;
	}
	return sum;
}
#endif	/* !USE_BLAS */

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
mmap_fd(int fd, size_t fz)
{
	void *p;

	p = mmap(NULL, fz, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
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
	struct stat st;
	glodfn_t res;

	if ((res.fd = open(fn, flags)) < 0) {
		;
	} else if (fstat(res.fd, &st) < 0) {
		res.fb = (glodf_t){.z = 0U, .d = NULL};
		goto clo;
	} else if ((res.fb = mmap_fd(res.fd, st.st_size)).d == NULL) {
	clo:
		close(res.fd);
		res.fd = -1;
	}
	return res;
}

static __attribute__((unused)) int
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

struct dl_file_s {
	uint8_t magic[4U];
	uint8_t flags[4U];
	size_t nvis;
	size_t nhid;

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
pump(const char *file)
{
	static struct dl_rbm_s res;
	static struct dl_rbm_priv_s p[1];

	if (UNLIKELY((p->f = mmap_fn(file, O_RDWR)).fd < 0)) {
		goto out;
	} else if (UNLIKELY(p->f.fb.z < sizeof(struct dl_file_s))) {
		goto out;
	}
	with (struct dl_file_s *fl = p->f.fb.d) {
		float_t *dp = fl->data + fl->off;

		res.nvis = fl->nvis;
		res.vbias = dp;
		dp += fl->nvis;

		res.nhid = fl->nhid;
		res.hbias = dp;
		dp += fl->nhid;

		res.w = dp;

		p->wtr = tr(res.w, fl->nvis, fl->nhid);
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

static dl_rbm_t
crea(const char *file, struct dl_file_s fs)
{
	static glodfn_t f;
	size_t z;
	size_t fz;
	int fd;

	if (UNLIKELY((fd = open(file, O_CREAT | O_RDWR | O_TRUNC, 0666)) < 0)) {
		return NULL;
	}
	/* compute total file size */
	z = fs.nvis + fs.nhid + fs.nvis * fs.nhid;
	ftruncate(fd, fz = z * sizeof(float) + sizeof(fs));

	if (UNLIKELY((f.fb = mmap_fd(f.fd = fd, fz)).d == NULL)) {
		goto out;
	}
	/* just copy the file header over */
	fs.off = 0U;
	memcpy(f.fb.d, &fs, sizeof(fs));

	munmap_fn(f);
	with (dl_rbm_t m = pump(file)) {
		const float vnois = .1f;
		const float hnois = .01f;
		const float wnois = 1.f / (m->nvis * m->nhid);

		/* wobble vbiasses */
		for (size_t i = 0; i < m->nvis; i++) {
			const float x = dr_rand_uni();
			m->vbias[i] = log(vnois * x);
		}

		/* wobble hbiasses */
		for (size_t j = 0; j < m->nhid; j++) {
			const float x = dr_rand_norm();
			m->hbias[j] = hnois * x;
		}

		/* wobble weights */
		for (size_t k = 0; k < m->nvis * m->nhid; k++) {
			const float x = dr_rand_norm();
			m->w[k] = wnois * x;
		}

		/* and wobble the transpose too */
		with (struct dl_rbm_priv_s *p = m->priv) {
			p->wtr = tr(m->w, m->nvis, m->nhid);
		}
		return m;
	}

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
			spsv = realloc(spsv, nu *= sizeof(*spsv));
			spsz = nu;
		}
		/* assign index/value pair */
		spsv[i].i = v;
		spsv[i].v = c;
		i++;
	}
	return (spsv_t){.z = i, .v = spsv};
}

static __attribute__((unused)) void
popul_ui8(float *restrict x, const uint8_t *n, size_t z)
{
	for (size_t i = 0; i < z; i++) {
		x[i] = (float)(int)n[i];
	}
	return;
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

		res += c;
		x[i] = (float)(int)c;
	}
	return res;
}

static __attribute__((unused)) float
poiss_lambda_ui8(const uint8_t *n, size_t z)
{
	/* N is the number of total incidents (occurence times count) */
	long unsigned int N = 0UL;

	/* calc N and n */
	for (size_t i = 0; i < z; i++) {
		N += n[i];
	}
	return (float)N / (float)z;
}

static __attribute__((unused)) float
poiss_lambda_f(const float *v, size_t z)
{
	/* incident number */
	long unsigned int N = 0UL;

	for (size_t i = 0; i < z; i++) {
		N += (long int)v[i];
	}
	return log((float)N) / (float)z;
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
static size_t N;

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
		h[j] = b[j] + cblas_sdot(nvis, w(j), 1U, vis, 1U);
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
		v[i] = b[i] + cblas_sdot(nhid, w(i), 1U, hid, 1U);
	}
#undef w
	return 0;
}

static ni int
expt_vis(float *restrict v, dl_rbm_t m, const float vis[static m->nvis])
{
	const size_t nvis = m->nvis;
#if defined SALAKHUTDINOV
	float nor = 0.f;
#endif	/* SALAKHUTDINOV */

	DEBUG(dump_layer("Va", vis, nvis));

#if defined SALAKHUTDINOV
	/* calc \sum exp(v) */
	for (size_t i = 0; i < nvis; i++) {
		nor += v[i] = exp(vis[i]);
	}
	with (const float norm = (float)N / nor) {
		for (size_t i = 0; i < nvis; i++) {
			v[i] = v[i] * norm;
		}
	}
#elif defined GEHLER
	for (size_t i = 0; i < nvis; i++) {
		v[i] = exp(vis[i]);
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
		v[i] = dr_rand_poiss(vis[i]);
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

static ni void
update_w(drbctx_t ctx)
{
	dl_rbm_t m = ctx->m;
	const float *vo = ctx->vo;
	const float *ho = ctx->ho;
	const float *vr = ctx->vr;
	const float *hr = ctx->hr;
	const size_t nv = m->nvis;
	const size_t nh = m->nhid;
#if !defined NDEBUG
	float mind = INFINITY;
	float maxd = -INFINITY;
#endif	/* !NDEBUG */

#define w(i, j)		m->w[i * nh + j]
#define wtr(i, j)	((struct dl_rbm_priv_s*)m->priv)->wtr[i * nv + j]
#define dw(i, j)	ctx->dw[i * nh + j]

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
	static void __upd(
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
	N = popul_sv(vo, nv, sv);

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
dream(drbctx_t ctx, spsv_t sv)
{
#define m	ctx->m
#define vo	ctx->vo
#define ho	ctx->ho
#define vr	ctx->vr
	const size_t nv = m->nvis;

	/* populate from input */
	N = popul_sv(vo, nv, sv);

	/* vhv gibbs */
	prop_up(ho, m, vo);
	expt_hid(ho, m, ho);
	smpl_hid(ho, m, ho);
	/* hv gibbs */
	prop_down(vr, m, ho);
	expt_vis(vr, m, vr);
	smpl_vis(vr, m, vr);

	for (size_t i = 0; i < nv; i++) {
		uint8_t vi = (uint8_t)(int)vr[i];

		if (UNLIKELY(vi)) {
			printf("%zu\t%u\n", i, (unsigned int)vi);
		}
	}
	return;
#undef m
#undef vo
#undef ho
#undef vr
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
#include "dbn-test.xh"
#include "dbn-test.x"
#if defined __INTEL_COMPILER
# pragma warning (default:593)
# pragma warning (default:181)
#endif	/* __INTEL_COMPILER */

int
main(int argc, char *argv[])
{
	struct glod_args_info argi[1];
	dl_rbm_t m = NULL;
	int res = 0;

	if (glod_parser(argc, argv, argi)) {
		res = 1;
		goto out;
	}

	init_rand();
	if (argi->create_given) {
		struct dl_file_s ini = {
			.nvis = 32768U + 4096U,
			.nhid = 256U,
		};

		/* just create the machine */
		if ((m = crea("test.rbm", ini)) == NULL) {
			res = 1;
		}
		goto wrout;
	}
	/* all other options need the machine file */
	/* read the machine file */
	if (UNLIKELY((m = pump("test.rbm")) == NULL)) {
		res = 1;
		goto wrout;
	}

	if (argi->check_given) {
		res = check(m);
		goto wrout;
	}

	/* from now on we're actually doing something with the machine */
	static struct drbctx_s ctx[1];
	const int fd = STDIN_FILENO;

	init_drbctx(ctx, m);
	if (argi->train_given) {
		for (spsv_t sv; (sv = read_tf(fd)).z; train(ctx, sv));
		final_update_w(ctx);
		final_update_b(ctx);
	} else if (argi->dream_given) {
		for (spsv_t sv; (sv = read_tf(fd)).z; dream(ctx, sv));
	}

	/* just to deinitialise resources */
	(void)read_tf(-1);
	fini_drbctx(ctx);
wrout:
	deinit_rand();
	dump(m);
out:
	glod_parser_free(argi);
	return res;
}
