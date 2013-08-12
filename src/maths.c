/*** maths.h -- calculation outsourcing
 *
 * Copyright (C) 2009-2013 Sebastian Freundt
 *
 * Author:  Sebastian Freundt <freundt@ga-group.nl>
 *
 * This file is part of drbang.
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
#if defined HAVE_CONFIG_H
# include "config.h"
#endif	/* HAVE_CONFIG_H */
#include <stdint.h>
#include <tgmath.h>
#include "maths.h"
#include "nifty.h"

/* tg defs */
#undef factorial
#undef poiss
#undef sigma
#undef softmax

#define PREFER_NUMERICAL_STABILITY_OVER_SPEED


float
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

double
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

long double
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

float
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

double
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

long double
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

float
sigmaf(float x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1.f + tanh(x / 2.f)) / 2.f;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1.f / (1.f + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

double
sigma(double x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1. + tanh(x / 2.)) / 2.;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1. / (1. + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

long double
sigmal(long double x)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	return (1.L + tanh(x / 2.L)) / 2.L;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	return 1.L / (1.L + exp(-x));
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

void
softmaxf(float *restrict tgt, const float *src, size_t z)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	float max = -INFINITY;
	float sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		if (src[i] > max) {
			max = src[i];
		}
	}
	for (size_t i = 0; i < z; i++) {
		sm += exp(src[i] - max);
	}
	with (const float lgsm = log(sm) + max) {
		for (size_t i = 0; i < z; i++) {
			tgt[i] = exp(src[i] - lgsm);
		}
	}
	return;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	float sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		sm += tgt[i] = exp(src[i]);
	}
	for (size_t i = 0; i < z; i++) {
		tgt[i] = tgt[i] / sm;
	}
	return;
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

void
softmax(double *restrict tgt, const double *src, size_t z)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	double max = -INFINITY;
	double sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		if (src[i] > max) {
			max = src[i];
		}
	}
	for (size_t i = 0; i < z; i++) {
		sm += exp(src[i] - max);
	}
	with (const double lgsm = log(sm) + max) {
		for (size_t i = 0; i < z; i++) {
			tgt[i] = exp(src[i] - lgsm);
		}
	}
	return;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	double sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		sm += tgt[i] = exp(src[i]);
	}
	for (size_t i = 0; i < z; i++) {
		tgt[i] = tgt[i] / sm;
	}
	return;
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

void
softmaxl(long double *restrict tgt, const long double *src, size_t z)
{
#if defined PREFER_NUMERICAL_STABILITY_OVER_SPEED
	long double max = -INFINITY;
	long double sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		if (src[i] > max) {
			max = src[i];
		}
	}
	for (size_t i = 0; i < z; i++) {
		sm += exp(src[i] - max);
	}
	with (const long double lgsm = log(sm) + max) {
		for (size_t i = 0; i < z; i++) {
			tgt[i] = exp(src[i] - lgsm);
		}
	}
	return;
#else  /* !PREFER_NUMERICAL_STABILITY_OVER_SPEED */
	long double sm = 0.f;

	for (size_t i = 0; i < z; i++) {
		sm += tgt[i] = exp(src[i]);
	}
	for (size_t i = 0; i < z; i++) {
		tgt[i] = tgt[i] / sm;
	}
	return;
#endif	/* PREFER_NUMERICAL_STABILITY_OVER_SPEED */
}

/* maths.c ends here */
