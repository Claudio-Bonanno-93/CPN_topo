#ifndef CPN_OP_H
#define CPN_OP_H

#include "macro.h"
#include "rng.h"
#include "endianness.h"

#include <complex.h>
#include <math.h>

typedef double complex cmplx; // short-hand for complex number with double precision real and imaginary part

// complex vector operations

// v=0
inline void vector_zero(cmplx *v)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif
	
	int i;
	for (i=0; i<N; i++)
	{
		v[i]=0.0+I*0.0;
	}
}

// v=w
inline void vector_equal(cmplx *v, cmplx const * const w)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	int i;
	for (i=0; i<N; i++)
	{
		v[i]=w[i];
	}
}

// v +=w
inline void vector_sum(cmplx *v, cmplx const * const w)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	__assume_aligned(&(w), DOUBLE_ALIGN);
	#endif

	int i;
	for (i=0; i<N; i++)
	{
		v[i]+=w[i];
	}
}

// res = a*v + b*w with a,b complex
inline void vector_linear_combination_cmplx_coeff(cmplx *res, cmplx const * const v, cmplx const * const w, cmplx const a, cmplx const b)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(res), DOUBLE_ALIGN);
	__assume_aligned(&(v), DOUBLE_ALIGN);
	__assume_aligned(&(w), DOUBLE_ALIGN);
	#endif

	int i;
	for (i=0; i<N; i++)
	{
		res[i] = a*v[i] + b*w[i];
	}
}

// v = a*v + b*w with a,b real
inline void vector_linear_combination_real_coeff(cmplx *v, cmplx const * const w, double const a, double const b)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	__assume_aligned(&(w), DOUBLE_ALIGN);
	#endif

	int i;
	for (i=0; i<N; i++)
	{
		v[i] = a*v[i] + b*w[i];
	}
}

// res = (v,w) = conj(v) \dot w = sum_{i=1}^{N} conj(v)[i] w[i]
inline cmplx vector_scalar_product(cmplx const * const v, cmplx const * const w)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	__assume_aligned(&(w), DOUBLE_ALIGN);
	#endif

	int i;
	cmplx res = 0.0 + I * 0.0;
	for( i=0 ; i<N ; i++)
	{
		res += conj(v[i]) * w[i];
	}

  return res;
}

// return |v|^2 = sum_{i=1}^N |v[i]|^2
inline double vector_norm(cmplx const * const v)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	return ( creal(vector_scalar_product(v,v)) );
}

// return |v|= sqrt{ sum_{i=1}^N |v[i]|^2 }
inline double vector_abs(cmplx const * const v)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	return ( sqrt(vector_norm(v)) );
}

// v *= constant
inline void vector_times_real_const(cmplx *v, double const constant)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	int i;
	for (i=0; i<N; i++)
	{
		v[i] *= constant;
	}
}

// v -> v/|v|
inline void vector_normalization(cmplx *v)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	double normalization=vector_abs(v);
	normalization=1.0/normalization;
	vector_times_real_const(v, normalization);
}

// complex number operations

// return |z|^2 = Re(z)^2 + Im(z)^2
inline double cmplx_norm(cmplx const z) 
{ 
	return ( creal(z)*creal(z) + cimag(z)*cimag(z) ); 
}

// return |z| = sqrt{ Re(z)^2 + Im(z)^2 }
inline double cmplx_abs(cmplx const z)
{
	return ( sqrt(cmplx_norm(z)) );
}

// return principal value of arg(z) = Im log(z) = atan2(Im(z),Re(z))
// if z = r exp(i theta), this function returns arg(z) = theta restricted the the interval [-pi, pi]
inline double arg(cmplx const z)
{
	double x, y;
	x=creal(z); 
	y=cimag(z); 
	return atan2(y,x);
}

// complex random

inline cmplx cmplx_rand_num(RNG_Param *rng_state)
{
	return( rand_num(rng_state) + I * rand_num(rng_state) );
}

inline void vector_random_cold(cmplx *v, RNG_Param *rng_state)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	int k;
	v[0] = 1.0 + I * 0.0;
	for(k=1; k<N; k++)
	{	
		// (2*r-1) * 10^(-6) = uniform random number in ( -10^(-6), 10^(-6) ) being r uniform random in (0,1)
		v[k] = ( (2.0 * cmplx_rand_num(rng_state) ) - (1.0 + I * 1.0) ) * (1.0e-06);
	}
	vector_normalization(v);
}

inline void vector_random_hot(cmplx *v, RNG_Param *rng_state)
{
	#ifdef __INTEL_COMPILER
	__assume_aligned(&(v), DOUBLE_ALIGN);
	#endif

	int k;
	for(k=0; k<N; k++)
	{
		v[k] = cmplx_rand_num(rng_state);
	}
	vector_normalization(v);
}

// read/write functions defined in lib/cpn_cmplx_op.c
void write_cmplx_num_on_file_binary_big_endian(FILE *, cmplx const * const);
void write_cmplx_num_on_file_binary_little_endian(FILE *, cmplx const * const);
void write_cmplx_num_on_file_binary(FILE *, cmplx const * const);
void write_cmplx_vec_on_file_binary(FILE *, cmplx const * const);

void read_cmplx_num_from_file_binary_big_endian(FILE *, cmplx *);
void read_cmplx_num_from_file_binary_little_endian(FILE *, cmplx *);
void read_cmplx_num_from_file_binary(FILE *, cmplx *);
void read_cmplx_vec_from_file_binary(FILE *, cmplx *);

#endif
