#ifndef ENDIANNESS_H
#define ENDIANNESS_H

#include "macro.h"
#include <stdlib.h>

// in lib/endianness.c
int is_little_endian(void);
void bytes_swap(void *p);

#endif
