#ifndef GEOM_H
#define GEOM_H

#include "cpn_param.h"

// geometry struct, stores near-neighbors relations
typedef struct Geometry{

	long ** up; // up[i][j] forward  near-neighbor of site i along direction j
	long ** dn; // dn[i][j] backward near-neighbor of site i along direction j

} Geometry;

// rectangle struct, stores rectangle sites
typedef struct Rectangle {

	long *rect_sites;
	int d_size_rect[2];
	long d_vol_rect;

} Rectangle;

// in lib/geometry.c
void init_geometry(Geometry *, CPN_Param const * const);
void si_to_cart(long *, long, CPN_Param const * const);
long cart_to_si(long const * const, CPN_Param const * const);
double square_distance(long const, long const, CPN_Param const * const param);
long periodic_condition(int const, int const);
long cart_to_si_rectangle(int const * const, Rectangle const * const);
void init_rectangle(Rectangle *, int const, CPN_Param const * const);
void init_rectangles_hierarchic_upd(Rectangle **, CPN_Param const * const);
void free_geometry(Geometry *, CPN_Param const * const);
void free_rectangles_hierarchic_upd(Rectangle *, CPN_Param const * const);
void free_rectangle(Rectangle *);

#endif
