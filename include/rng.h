#ifndef RNG_H
#define RNG_H

#include "cpn_param.h"

#include <math.h>
#include <float.h> // to define macro DBL_EPSILON ~ 2.2 x 10^(-16)

// ran2 constants
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

// ran 2 rng state struct
typedef struct RNG_Param {

	long idum;
	long idum2;
	long iy;
	long iv[NTAB];

} RNG_Param;

// in lib/rng.c

// ran2 functions
void ran2_std_init_rng_state(RNG_Param *, CPN_Param const * const);
double ran2_rng(RNG_Param *);
void init_ran2_rng_state(RNG_Param *, CPN_Param const * const);
void read_ran2_rng_state_from_file(RNG_Param *, CPN_Param const * const);
void write_ran2_rng_state(FILE *, RNG_Param const * const);

// general random functions

// polar angle RNG for over-heat-bath update of CP^(N-1) Models
double rng_theta_U(RNG_Param *, double const);
double rng_theta_z(RNG_Param *, double const);

// use this function to call a uniform Random Number Generator (RNG) in (0,1) with your favorite algorithm
double rand_num(RNG_Param *);

// use these functions to read, write and initialize rng state
void init_rng_state(RNG_Param *, CPN_Param const * const);
void write_rng_state(RNG_Param const * const, CPN_Param const * const);

#endif
