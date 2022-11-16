#ifndef GEO_C
#define GEO_C

#include "../include/geometry.h"

// create geometry, defining up and down near-neighbors
void init_geometry(Geometry * geo, CPN_Param const * const param)
{
	int mu, err;
	long cart_coord[2];
	long i, value, aux;

	// allocate memory
	err=posix_memalign((void**)&(geo->up), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	err=posix_memalign((void**)&(geo->dn), (size_t)INT_ALIGN, (size_t) param->d_volume * sizeof(long *));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
	exit(EXIT_FAILURE);
	}

	for(i=0; i<(param->d_volume); i++)
	{
		err=posix_memalign((void**)&(geo->up[i]), (size_t)INT_ALIGN, (size_t) 2 * sizeof(long));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		err=posix_memalign((void**)&(geo->dn[i]), (size_t)INT_ALIGN, (size_t) 2 * sizeof(long));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the geometry! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	// initialize
	for (i=0; i<param->d_volume; i++) 
	{
		si_to_cart(cart_coord,i,param);
		for(mu=0; mu<2; mu++)
		{
			value=cart_coord[mu]; // site component along mu direction
			
			// up near-neighbor
			aux=value+1;
			if(aux >= param->d_size[mu]) aux-=param->d_size[mu];
			cart_coord[mu]=aux;
			aux=cart_to_si(cart_coord, param);
			geo->up[i][mu]=aux;

			// down near-neighbor
			aux=value-1;
			if(value<0) aux+=param->d_size[i];
			cart_coord[mu]=aux;
			aux=cart_to_si(cart_coord, param);
			geo->dn[i][mu]=aux;

			cart_coord[mu]=value;
		} 
	}
}	

// convert from single index (lexicographic) to cartesian coordinates
void si_to_cart(long *cart_coord, long r, CPN_Param const * const param)
{
	cart_coord[0]=(r+param->d_size[0])%param->d_size[0];
	cart_coord[1]=((r/param->d_size[0])+param->d_size[1])%param->d_size[1];
}

// convert from cartesian to single index (lexicographic) coordinates
long cart_to_si(long const * const cart_coord, CPN_Param const * const param)
{
	long r;
	r=(cart_coord[0]+param->d_size[0])%param->d_size[0] + ((cart_coord[1]+param->d_size[1])%param->d_size[1])*param->d_size[0];
	return r;
}

// assume boundary conditions on [0, L_max - 1] and compute the value of coord
long periodic_condition(int const coord, int const L_max)
{
	if(coord<0) return (L_max+coord%L_max);
	else if(coord>=L_max) return (coord%L_max);
	else return coord;
}

// compute the square distance between points i and j on a periodic lattice (toroidal geometry)
double square_distance(long const i, long const j, CPN_Param const * const param)
{
	int mu;
	long x[2], y[2];
	double d[2], half_size[2];
	double distance2=0.0;

	si_to_cart(x, i, param); // i --> x
	si_to_cart(y, j, param); // j --> y
	for (mu=0; mu<2; mu++)
	{
		d[mu]=( (double) (labs(x[mu] - y[mu])) ); // d_mu = | x_mu - y_mu |
		half_size[mu] = ((double) param->d_size[mu])/2.0;
		if ( d[mu] > half_size[mu] ) d[mu] = ((double) param->d_size[mu]) - d[mu]; // periodic boundaries: if d_mu > L/2 ==> d_mu = L - d_mu
		distance2 += d[mu]*d[mu]; // distance^2 = sum_mu d_mu^2
	}
	return distance2;
}

// convert from cartesian to single index (lexicographic) coordinates for a given rectangle
long cart_to_si_rectangle(int const * const cart_coord, Rectangle const * const most_update)
{	
	int mu;
	long ris=0, aux=1;

	for(mu=0; mu<2; mu++)
	{
		ris+=cart_coord[mu]*aux;
		aux*=most_update->d_size_rect[mu];
	}
	return ris;
}

void init_rectangle(Rectangle *most_update, int const L_R, CPN_Param const * const param)
{ 
	// compute dimensions of rectangle and ranges of rect coordinates
	int aux_L[2], size_min[2], size_max[2];
	int mu,err;

	aux_L[0]=2*L_R+1;
	aux_L[1]=2*L_R+param->d_L_defect;

	size_min[0]=param->d_size[0]-1-L_R;
	size_min[1]=-L_R;

	size_max[0]=param->d_size[0]+L_R;
	size_max[1]=param->d_L_defect+L_R;

	// d_size_rect[i] must not exceed d_size[i]
	// if so, the exceeding dimension of the rectangle is just the respective dimension of the whole lattice
	// and the i-th coordinate just ranges from 0 to d_size[i]
	for (mu=0; mu<2; mu++)
	{
		if(aux_L[mu]>=param->d_size[mu])
		{
			aux_L[mu]=param->d_size[mu];
			size_min[mu]=0;
			size_max[mu]=param->d_size[mu];
		}
	}	

	// store the volume of the rectangle
	long V_rect=1;
	for(mu=0; mu<2; mu++) V_rect*=aux_L[mu];

	// allocate rectangle
	err=posix_memalign((void **) &(most_update->rect_sites), (size_t) INT_ALIGN, (size_t) V_rect * sizeof(long));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating a rectangle!\n");
		exit(EXIT_FAILURE);
	}

	// save dimensions and volume of the rectangle
	for(mu=0; mu<2; mu++) most_update->d_size_rect[mu]=aux_L[mu];
	most_update->d_vol_rect=V_rect;

	// save lexicographical index of sites of the rectangle
	int coord[2];		// cartesian coordinates on the whole lattice
	long real_coord[2];	// cartesian coordinates after using periodic conditions
	long r,r_rect;		// lexicographical index on the whole lattice and lexicographical index on the rectangle
	int coord_rect[2];	// cartesian coordinates on the rectangle

	// loop on the rectangle coordinates
	coord_rect[0]=0;
	for(coord[0]=size_min[0];coord[0]<size_max[0];coord[0]++)
	{
		coord_rect[1]=0;
		for(coord[1]=size_min[1];coord[1]<size_max[1];coord[1]++)
		{
			// compute the real coordinates on the periodic lattice
			for(mu=0; mu<2; mu++) real_coord[mu]=periodic_condition(coord[mu],param->d_size[mu]);

			r_rect=cart_to_si_rectangle(coord_rect,most_update);	// from cartesian coordinates on the rectangle to single index on the rectangle
			r=cart_to_si(real_coord, param);	// from cartesian coordinates on the lattice to single index on the lattice
			most_update->rect_sites[r_rect]=r; // map the single index coordinate on the rectangle to its corresponding value on the whole lattice
			
			coord_rect[1]++;
		} // end of loop for mu=1 direction
		coord_rect[0]++;
	} // end of loop for mu=0 direction	
}
  
void init_rectangles_hierarchic_upd(Rectangle **most_update, CPN_Param const * const param)
{
	if ( param->d_N_hierarc_levels == 0 ) (void) most_update; // to suppress compiler warning about unused variable
	else
	{
		int i,err;
		err=posix_memalign( (void **) most_update, (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(Rectangle) );
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating rectangles!\n");
			exit(EXIT_FAILURE);
		}
		for(i=0; i<param->d_N_hierarc_levels; i++) init_rectangle(&((*most_update)[i]), param->d_L_rect[i],param);
	}
}
 
// free memory allocated for geometry
void free_geometry(Geometry *geo, CPN_Param const * const param)
{
	long r;
	for(r=0; r<param->d_volume; r++)
	{
		free(geo->up[r]);
		free(geo->dn[r]);
	}
	free(geo->up);
	free(geo->dn);
}

// free memory allocated for rectangles
void free_rectangles_hierarchic_upd(Rectangle *most_update, CPN_Param const * const param)
{
	if ( param->d_N_hierarc_levels == 0 ) (void) most_update; // to suppress compiler warning about unused variable
	else
	{
		int i;
		for(i=0;i<param->d_N_hierarc_levels;i++) free_rectangle(&(most_update[i]));
		free(most_update);
	}
}

void free_rectangle(Rectangle *most_update)
{
	free(most_update->rect_sites);
}

#endif
