#ifndef ENDIANNESS_C
#define ENDIANNESS_C

#include "../include/endianness.h"

// check if machine supports little or big endian. 0=little endian, 1=big endian
int is_little_endian(void)
{
	int i = 1;
	char *c = ( (char *) (&i) );

	if (c[0] == 1) return 0; // little endian
	else return 1; // big endian
}

// swap bytes for double variable
void bytes_swap(void *p)
{
	char *c = p;
	size_t low, high;
	char temp;

	for(low=0, high=sizeof(double)-1; high>low; low++, high--)
	{
		temp=c[low];
		c[low] = c[high];
		c[high] = temp;
	}
}

#endif
