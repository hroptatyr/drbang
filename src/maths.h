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
#if !defined INCLUDED_maths_h_
#define INCLUDED_maths_h_

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <tgmath.h>

extern float factorialf(uint8_t n);
extern double factorial(uint8_t n);
extern long double factoriall(uint8_t n);

extern float poissf(float lambda, uint8_t n);
extern double poiss(double lambda, uint8_t n);
extern long double poissl(long double lambda, uint8_t n);
#define poiss(x, n)	__TGMATH_BINARY_FIRST_REAL_ONLY(x, n, poiss)

extern float sigmaf(float x);
extern double sigma(double x);
extern long double sigmal(long double x);
#define sigma(x)	__TGMATH_UNARY_REAL_ONLY(x, sigma)

extern void softmaxf(float *restrict tgt, const float *src, size_t z);
extern void softmax(double *restrict tgt, const double *src, size_t z);
extern void softmaxl(long double *restrict tgt, const long double *src, size_t);
#define softmax(tgt, src, z)	__TGMATH_ARRAY_REAL_ONLY(tgt, src, z, softmax)

#define __TGMATH_ARRAY_REAL_ONLY(tgt, src, z, fct)			\
	(__extension__((sizeof(*src) == sizeof(double) ||		\
			__builtin_classify_type(*src) != 8)		\
		       ? (void)fct((void*)tgt, (const void*)src, z)	\
		       : (sizeof(*src) == sizeof(float))		\
		       ? (void)fct##f((void*)tgt, (const void*)src, z)	\
		       : (void)__tgml(fct)((void*)tgt, (const void*)src, z)))

#endif	/* INCLUDED_maths_h_ */
