#ifndef CPN_OP_C
#define CPN_OP_C

#include "../include/cpn_cmplx_op.h"

// inline functions defined in include/cpn_cmplx_op.h
void vector_zero(cmplx *v);
void vector_equal(cmplx *v, cmplx const * const w);
void vector_sum(cmplx *v, cmplx const * const w);
void vector_linear_combination_cmplx_coeff(cmplx *res, cmplx const * const v, cmplx const * const w, cmplx const a, cmplx const b);
void vector_linear_combination_real_coeff(cmplx *v, cmplx const * const w, double const a, double const b);
cmplx vector_scalar_product(cmplx const * const v, cmplx const * const w);
double vector_norm(cmplx const * const v);
double vector_abs(cmplx const * const v);
void vector_times_real_const(cmplx *v, double const constant);
void vector_normalization(cmplx *v);
double cmplx_norm(cmplx const z);
double cmplx_abs(cmplx const z);
double arg(cmplx const z);
cmplx cmplx_rand_num(RNG_Param *rng_state);
void vector_random_cold(cmplx *v, RNG_Param *rng_state);
void vector_random_hot(cmplx *v, RNG_Param *rng_state);

// read/write functions

// write
void write_cmplx_num_on_file_binary_big_endian(FILE *fp, cmplx const * const a)
{
	size_t err=0;
	double aux;

	// write real part of a
	aux=creal(*a);
	err+=fwrite(&(aux), sizeof(double), 1, fp);
	// write imag part of a
	aux=cimag(*a);
	err+=fwrite(&(aux), sizeof(double), 1, fp);

	if(err!=2)
	{
		fprintf(stderr, "Problem in writing a complex number on binary file (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
}

void write_cmplx_num_on_file_binary_little_endian(FILE *fp, cmplx const * const a)
{
	size_t err=0;
	double aux;

	// write real part of a
	aux=creal(*a);
	bytes_swap(&aux);
	err+=fwrite(&(aux), sizeof(double), 1, fp);
	// write real part of a
	aux=cimag(*a);
	bytes_swap(&aux);
	err+=fwrite(&(aux), sizeof(double), 1, fp);

	if(err!=2)
	{
		fprintf(stderr, "Problem in writing a complex number on binary file (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
}

void write_cmplx_num_on_file_binary(FILE *fp, cmplx const * const a)
{
	// if machine is little endian, write after bytes swap
	if (is_little_endian() == 0) write_cmplx_num_on_file_binary_little_endian(fp, a);
	// if machine is big endian, write without bytes swap
	else write_cmplx_num_on_file_binary_big_endian(fp, a);
}

void write_cmplx_vec_on_file_binary(FILE *fp, cmplx const * const v)
{
	int k;
	for (k=0; k<N; k++) write_cmplx_num_on_file_binary(fp, &(v[k]));
}


void read_cmplx_num_from_file_binary_big_endian(FILE *fp, cmplx *a)
{
	size_t err=0;
	double real, imag;

	// read real and imag parts
	err+=fread(&(real), sizeof(double), 1, fp);
	err+=fread(&(imag), sizeof(double), 1, fp);

	if(err!=2)
	{
		fprintf(stderr, "Problem in reading a complex number from binary file (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else (*a)=real+I*imag;
}

void read_cmplx_num_from_file_binary_little_endian(FILE *fp, cmplx *a)
{
	size_t err=0;
	double real, imag;

	// read real and imag parts
	err+=fread(&(real), sizeof(double), 1, fp);
	err+=fread(&(imag), sizeof(double), 1, fp);

	if(err!=2)
	{
		fprintf(stderr, "Problem in reading a complex number from binary file (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		bytes_swap( (void *) (&real) );
		bytes_swap( (void *) (&imag) );
		(*a)=real+I*imag;
	}
}

void read_cmplx_num_from_file_binary(FILE *fp, cmplx *a)
{
	// if machine is little endian, write after bytes swap
	if (is_little_endian() == 0) read_cmplx_num_from_file_binary_little_endian(fp, a);
	// if machine is big endian, write without bytes swap
	else read_cmplx_num_from_file_binary_big_endian(fp, a);
}

void read_cmplx_vec_from_file_binary(FILE *fp, cmplx *v)
{
	int k;
	for (k=0; k<N; k++) read_cmplx_num_from_file_binary(fp, &(v[k]));
}

#endif
