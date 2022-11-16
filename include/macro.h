#ifndef MACRO_H
#define MACRO_H

#include "../config.h"

// CP^(N-1) Models number of colors N
#define N 2

// smallest number
static const double NUM_EPS=1.0e-15;

// pi constant
static const double pi=3.141592653589793238462643383279502884197169399375105820974944;

// improvement coefficients appearing in the Symanzik action
static const double c1=(4.0/3.0); 
static const double c2=(-1.0/12.0);

// memory alignements
#define INT_ALIGN 16
#define DOUBLE_ALIGN 32

// max length of unknown string
#define STD_STRING_LENGTH 50

// to activate posix_memalign in stdlib.h
#define _POSIX_C_SOURCE 200809L

#endif
