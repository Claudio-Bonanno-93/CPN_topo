#ifndef RNG_C
#define RNG_C

#include "../include/rng.h"

// Over-heat-bath polar angles rng. These functions are written independently of the particular uniform rng used

// Random extraction of polar angle for heat-bath update of gauge field U_mu(x) using Von Neumann algorithm
// Taken from "Monte Carlo simulation of CP^{N-1} models", M. Campostrini, P. Rossi, E. Vicari, Phys. Rev. D 46 (1992) 2647-2662
double rng_theta_U(RNG_Param *rng_state, double const a)
{
	double theta_test;
	if(a==0)
	{
		theta_test = pi * rand_num(rng_state); // theta random uniform in (0,pi)
	}
	else
	{
		double c, k0, y, x_test, eta, test;
		c = sqrt(a/2.0); 
		if (a<0.8) eta=0.73; 
		else eta=0.99; 
		k0=atan(pi*c); 
		do // P_true(theta) = Z exp( a cos(theta) ) is sampled using Von Neumann algorithm
		{ 
			x_test = rand_num(rng_state); // uniform random in (0,1)
			// theta_test extracted as a Lorentz distribution P_Lorentz(theta_test) = (c/k0) 1/(1+(c*theta_test)^2) with theta in (0,pi)
			theta_test = tan(k0*x_test)/c;
			// test = P_true(theta_test) / ( M P_Lorentz(theta_test) ), where M = k0 * Z * exp(2a) / ( c eta )
			test = eta*( 1.0 + c*c*theta_test*theta_test )*exp( ( cos(theta_test) - 1.0 ) * a ); 
			y = rand_num(rng_state); 
		} 
		while (y > test); // Von-Neumann test
	}
	return theta_test;
}

// Random extraction of polar angle for heat-bath update of scalar field z(x) using Von Neumann algorithm
// Taken from "Monte Carlo simulation of CP^{N-1} models", M. Campostrini, P. Rossi, E. Vicari, Phys. Rev. D 46 (1992) 2647-2662
double rng_theta_z(RNG_Param *rng_state, double const a)
{
	double R, exponent, cos_theta_max, sin_theta_max, theta_max, c, k1, k2, b0;
	double y, x_test, theta_test, test, eta;
	eta=0.99; 

	R = ((double)(N-1))/a; 
	exponent = (double)((2*N)-2);

	cos_theta_max = sqrt(1.0+R*R) - R ; // sin(theta_max)
	sin_theta_max = sqrt(1.0-cos_theta_max*cos_theta_max); // cos(theta_max)
  
	theta_max = acos(cos_theta_max); // theta_max = arccos( cos(theta_max) )
	c = sqrt( a*(cos_theta_max + R) ); 
	k1 = atan(c*(pi-theta_max)); 
	k2 = atan(c*theta_max); 

	do // P_true(theta) = Z exp( a cos(theta) ) sin^(2N-2)(theta) is sampled using Von Neumann algorithm
	{
		x_test = rand_num(rng_state); // uniform random in (0,1)
		theta_test = theta_max + tan( x_test*k1 + (x_test-1.0)*k2 )/c; // theta_test extracted as a Lorentz distribution centered around theta_max 
		b0 = sin(theta_test)/sin_theta_max;
		// test = P_true(theta_test) / (M P_lorentz(theta_test) ), where M = (k1 + k2) Z exp( 2a cos(theta_max) ) sin^(4N-4)(theta_max) / (c eta)
		test = eta*( pow(b0,exponent) )*( 1.0 + c*c*(theta_test-theta_max)*(theta_test-theta_max) )*exp( a * ( cos(theta_test) - cos_theta_max ) );
		y = rand_num(rng_state);
	}
	while (y>test); // Von-Neumann test
  
  return theta_test;
}

// When changing rng algorithm, just change these functions in this file

// Extract a uniform random number in (0,1) without endpoints
double rand_num(RNG_Param *rng_state)
{
	return (ran2_rng(rng_state));
}

// initialize rng state
void init_rng_state(RNG_Param *rng_state, CPN_Param const * const param)
{
	init_ran2_rng_state(rng_state, param);
}

// print rng state
void write_rng_state(RNG_Param const * const rng_state, CPN_Param const * const param)
{
	FILE *fp;

	fp=fopen(param->d_rng_file, "w"); // open rng state file
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_rng_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		// print seed
		fprintf(fp, "%ld\n", param->d_seed);
		// print rng state
		write_ran2_rng_state(fp, rng_state);
	}
	fclose(fp);
}

//--------------- ran2 rng functions --------------------

// Uniform random number generator "ran2" adapted from "Numerical Recipes in C - The Art of Scientific Computing (Second Edition)"
// W. H. Press, S. A.  Teukolsky, W. T. Vetterling and B. P. Flannery, Cambridge University Press, 1992
// ran2 returns a random number x uniformely distributed in (0,1) (without endpoints)
// idum must be initialized with a NON-ZERO NEGATIVE INTEGER SEED
// This version is modified to separate the initialization of the rng state from the generation of the random number
// Moreover, it is modified to return a double as in, e.g., https://github.com/PrincetonUniversity/athena-public-version/blob/master/src/utils/ran2.cpp

// standard initalization of ran2 rng state (see Numerical Recipes in C)
void ran2_std_init_rng_state(RNG_Param *rng_state, CPN_Param const * const param)
{
	int i, j;
	long k;

	rng_state->idum=param->d_seed;
	if (rng_state->idum > 0) (rng_state->idum) *= -1; // if seed is positive, make it negative
	rng_state->idum2=123456789;
	rng_state->iy=0;
	for (i=0; i<NTAB; i++) rng_state->iv[i]=0;

	if (rng_state->idum <= 0) 
	{ 
		if (-(rng_state->idum) < 1) rng_state->idum=1; 
		else rng_state->idum = -(rng_state->idum); 
		rng_state->idum2=(rng_state->idum); 
		for (j=NTAB+7;j>=0;j--) 
		{ 
			k=(rng_state->idum)/IQ1; 
			rng_state->idum=IA1*(rng_state->idum-k*IQ1)-k*IR1; 
			if (rng_state->idum < 0) rng_state->idum += IM1; 
			if (j < NTAB) rng_state->iv[j] = rng_state->idum; 
		} 
		rng_state->iy=rng_state->iv[0]; 
	}
}

// ran2 rng
double ran2_rng(RNG_Param *rng_state)
{
	int j;
	long k;
	double temp;

	k=(rng_state->idum)/IQ1; 
	rng_state->idum=IA1*(rng_state->idum-k*IQ1)-k*IR1; 
	if (rng_state->idum < 0) rng_state->idum += IM1; 
	k=rng_state->idum2/IQ2; 
	rng_state->idum2=IA2*(rng_state->idum2-k*IQ2)-k*IR2; 
	if (rng_state->idum2 < 0) rng_state->idum2 += IM2; 
	j=(int)(rng_state->iy/NDIV); 
	rng_state->iy=rng_state->iv[j]-rng_state->idum2; 
	rng_state->iv[j] = rng_state->idum; 
	if (rng_state->iy < 1) rng_state->iy += IMM1;

	if ( ( temp = AM * ((double)rng_state->iy) ) > RNMX) return RNMX;
	else return temp;
}

// initialize ran2 rng state
void init_ran2_rng_state(RNG_Param *rng_state, CPN_Param const * const param)
{
	if (param->d_rng_start == 0 ) ran2_std_init_rng_state(rng_state, param); // standard initialization from scratch
	if (param->d_rng_start == 1 ) read_ran2_rng_state_from_file(rng_state, param); // read last ran2 rng state
}

// read ran2 rng state from file
void read_ran2_rng_state_from_file(RNG_Param *rng_state, CPN_Param const * const param)
{
	long stored_seed, stored_idum, stored_idum2, stored_iy, stored_iv;
	int i, err;
	FILE * fp;

	// open rng state file for reading
	fp=fopen(param->d_rng_file, "r"); // open rng state file
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_rng_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	err = fscanf(fp, "%ld %ld %ld %ld", &stored_seed, &stored_idum, &stored_idum2, &stored_iy);
	if (err == 4)
	{
		if (stored_seed != param->d_seed)
		{
			fprintf(stderr, "Error: input file seed and rng state file seed are different (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		else
		{
			rng_state->idum=stored_idum;
			rng_state->idum2=stored_idum2;
			rng_state->iy=stored_iy;
		}
	}
	else
	{
		fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_rng_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	
	for (i=0; i<NTAB; i++)
	{
		err = fscanf(fp, "%ld", &stored_iv);
		if (err == 1) rng_state->iv[i]=stored_iv;
		else
		{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_rng_file, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}
}

// write ran2 rng state to file
void write_ran2_rng_state(FILE *fp, RNG_Param const * const rng_state)
{
	int i;
	fprintf(fp, "%ld %ld %ld ", rng_state->idum, rng_state->idum2, rng_state->iy);
	for (i=0; i<NTAB; i++) fprintf(fp, "%ld ", rng_state->iv[i]);
}

#endif
