#ifndef CPN_UPDATE_C
#define CPN_UPDATE_C

#include "../include/cpn_conf.h"

// perform a single step of Hasenbusch parallel tempering with hierarchic updates and translations
// Adapted from "Fighting topological freezing in the two-dimensional CP^{N-1} model", M. Hasenbusch, Phys. Rev. D 96, 054504 (2017), 1706.04443
void parallel_tempering_with_hierarchic_update(CPN_Conf * conf, Rectangle const * const most_update, Acc_Swap *swap_counter, CPN_Param const * const param,
												Geometry const * const geo, CPN_Conf *aux_conf, RNG_Param *rng_state)
{
	int i, begin_hierarc_level=0;
	for(i=0;i<param->d_N_replica_pt;i++)
	{
		/*------- Single parallel tempering step -------*/
		update_with_defect(&(conf[i]), geo, param, rng_state);
		if(param->d_N_replica_pt>1)
		{
			swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
			if ( param->d_N_hierarc_levels > 0 )
			{
				hierarchic_update_rectangle(conf, geo, param, i, begin_hierarc_level, most_update, swap_counter, aux_conf, rng_state);
			}
		}
		/*----------------------------------------------*/

		// increase counter of update steps
		conf[i].update_index++;
	}
}
	
// perform a hierarchic update on rectangles for conf[conf_label]
void hierarchic_update_rectangle(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, int const conf_label, int const hierarc_level,
                                 Rectangle const * const most_update, Acc_Swap *swap_counter, CPN_Conf *aux_conf, RNG_Param *rng_state)
{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
	{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
		{
			update_rectangle(&(conf[conf_label]), geo, param, &(most_update[hierarc_level]), rng_state);
			if(param->d_N_replica_pt>1) swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
		}
	} // end if
	else
	{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
		{
			update_rectangle(&(conf[conf_label]), geo, param, &(most_update[hierarc_level]), rng_state);   
			if(param->d_N_replica_pt>1) swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
			hierarchic_update_rectangle(conf, geo, param, conf_label, (hierarc_level+1), most_update, swap_counter, aux_conf, rng_state);
		}
	} // end else
}
	
// update only on a given rectangle in the presence of a defect (to be used during parallel tempering)
void update_rectangle(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, Rectangle const * const most_update, RNG_Param *rng_state)
{
	int j;
	for (j=0; j<param->d_num_micro; j++) microcanonic_sweep_rectangle(conf,geo,param,most_update); 
	overheatbath_sweep_rectangle(conf,geo,param,most_update,rng_state); 
}
	
// perform an update of the whole lattice in the presence of a defect
void update_with_defect(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, RNG_Param *rng_state)
{
	int j;
	for (j=0; j<param->d_num_micro; j++) microcanonic_sweep_lattice(conf,geo,param);
	overheatbath_sweep_lattice(conf,geo,param,rng_state);
}

// compute forward staple relative to link (i,mu)
cmplx staple_up(CPN_Conf const * const conf, Geometry const * const geo, long const i, int mu)
{
	int nu = 1 - mu;
	return (conf->U[geo->up[i][mu]][nu]*conj(conf->U[geo->up[i][nu]][mu])*conj(conf->U[i][nu]));
}

// compute backward staple relative to link (i,mu)
cmplx staple_down(CPN_Conf const * const conf, Geometry const * const geo, long const i, int mu)
{
	int nu = 1 - mu;
	return (conj( conf->U[geo->dn[geo->up[i][mu]][nu]][nu] )*conj(conf->U[geo->dn[i][nu]][mu])*conf->U[geo->dn[i][nu]][nu]);
}

// compute F_U = dS/dU relative to field U on link (i,k)
cmplx force_U(CPN_Conf const * const conf, Geometry const * const geo, CPN_Param const * const param, long const i, int const mu) 
{
	double sign;
	cmplx a, b, c, F_topo;

	// non-topological contributions to the force
	a=vector_scalar_product(conf->z[ geo->up[i][mu]], conf->z[i]);
	b=vector_scalar_product(conf->z[ geo->up[ geo->up[i][mu] ][mu]], conf->z[i]);
	c=vector_scalar_product(conf->z[ geo->up[i][mu] ], conf->z[ geo->dn[i][mu] ]);

	a *= conf->C[i][mu];
	b *= conf->C[i][mu] * conf->C[geo->up[i][mu]][mu] * conj(conf->U[geo->up[i][mu]][mu]);
	c *= conf->C[geo->dn[i][mu]][mu] * conf->C[i][mu] * conj(conf->U[geo->dn[i][mu]][mu]);

	// topological contribution to the force
	if (mu==0) {sign=1.0;}
	else {sign=-1.0;}

	F_topo = conj( staple_up(conf, geo, i, mu) - staple_down(conf, geo, i, mu) );
	F_topo *= I * sign * param->d_theta / (4.0 * pi * param->d_beta * ((double)N));

	return (c1*a +  c2*(b + c) + F_topo); // F_U
} 


// compute F_z = dS/dz relative to field z on site i
void force_z(CPN_Conf const * const conf, Geometry const * const geo, long i, cmplx *F_z) 
{
	int mu;
	cmplx coeff1, coeff2;
	cmplx x1[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx x2[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx aux[N] __attribute__((aligned(DOUBLE_ALIGN)));

	vector_zero(F_z); // F_z=0
	vector_zero(x1);  // x1=0
	vector_zero(x2);  // x2=0
	for(mu=0; mu<2; mu++)
	{	
		// contribution from first term of Symanzik-improved action

		// coeff1 = ( C*conj(U) )(i-mu)_mu
		coeff1=conf->C[geo->dn[i][mu]][mu] * conj(conf->U[geo->dn[i][mu]][mu]);
		// coeff2 = (C*U)(i)_mu
		coeff2=conf->C[i][mu] * conf->U[i][mu];
		// aux = coeff1 * z(i-mu) + coeff2 * x(i+mu)
		vector_linear_combination_cmplx_coeff(aux, conf->z[geo->dn[i][mu]], conf->z[geo->up[i][mu]], coeff1, coeff2);
		vector_sum(x1, aux); // x1 += aux

		// contribution from second term of Symanzik-improved action

		// coeff1 = (C*U)(i-2mu)_mu * ( C*conj(U) )(i-mu)_mu
		coeff1=conf->C[geo->dn[geo->dn[i][mu]][mu]][mu] * conf->C[geo->dn[i][mu]][mu] * conj(conf->U[geo->dn[i][mu]][mu] * conf->U[geo->dn[geo->dn[i][mu]][mu]][mu]);
		// coeff2 = (C*U)(i)_mu * (C*U)(i+mu)_mu
		coeff2=conf->C[geo->up[i][mu]][mu] * conf->C[i][mu] * conf->U[i][mu] * conf->U[geo->up[i][mu]][mu];
		// aux = coeff1  * z(i-2mu) + coeff2 * z(i+2mu)
		vector_linear_combination_cmplx_coeff(aux, conf->z[geo->dn[geo->dn[i][mu]][mu]], conf->z[geo->up[geo->up[i][mu]][mu]], coeff1, coeff2);
		vector_sum(x2, aux); // x2 +=aux
	}

	// compute total force
	vector_times_real_const(x1, c1); // x1 -> c1 x1
	vector_times_real_const(x2, c2); // x1 -> c2 x2
	vector_sum(F_z, x1); // F_z = c1 x1
	vector_sum(F_z, x2); // F_z = c1 x1 + c2 x2
}

void microcanonic_sweep_rectangle(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, Rectangle const * const most_update)
{
	long i, n;
	int mu;
	for (i=0; i<(most_update->d_vol_rect); i++)
	{
		n=most_update->rect_sites[i];
		for (mu=0; mu<2; mu++) microcanonic_single_link_U(conf, geo, param, n, mu);
		microcanonic_single_site_z(conf, geo, n);
	}
}	

void microcanonic_sweep_lattice(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param)
{
	long i;
	int mu;
	for (i=0; i<param->d_volume ; i++)
	{
		for (mu=0; mu<2; mu++) microcanonic_single_link_U(conf, geo, param, i, mu);
		microcanonic_single_site_z(conf, geo, i); 
	}
}

// over-relaxation update of U and z fields
// Taken from "Monte Carlo simulation of CP^{N-1} models", M. Campostrini, P. Rossi, E. Vicari, Phys. Rev. D 46 (1992) 2647-2662

// over-relaxation update of field U on link (i,mu)
void microcanonic_single_link_U(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, long const i, int const mu) 
{ 
	double p, mod_F_U; 
	cmplx F_U;

	F_U = force_U(conf, geo, param, i, mu);
	p = creal( conj(conf->U[i][mu]) * F_U ); // project U along F_U for link (i, mu)
	mod_F_U=cmplx_norm(F_U); // |U|^2

	conf->U[i][mu] *= -1.0; // U_new = -U_old if F_U=0
	if (mod_F_U!=0) conf->U[i][mu] += ( 2.0 * p * F_U / mod_F_U ); // U_new = 2 p U_old F_U/|F_U|^2 - U_old if F_U!=0
}

// over-relaxation of field z on site i
void microcanonic_single_site_z(CPN_Conf * conf, Geometry const * const geo, long const i) 
{
	double p, mod_F_z, b; 
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx a;

	force_z(conf, geo, i, F_z);
	a=vector_scalar_product(conf->z[i], F_z); // project z along F_z for site i
	p = creal(a);
	mod_F_z = vector_norm(F_z); // |F_z|^2
	b = (2.0 * p / mod_F_z);

	if (mod_F_z!=0) vector_linear_combination_real_coeff(conf->z[i], F_z, -1.0, b); // z_new = -z_old + b F_z if F_z!=0
	else vector_times_real_const(conf->z[i], -1.0); // z_new = -z_old if F_z=0
} 

void overheatbath_sweep_rectangle(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, Rectangle const * const most_update, RNG_Param *rng_state)
{
	long i, n;
	int mu;
	for (i=0; i<most_update->d_vol_rect; i++)
	{
		n=most_update->rect_sites[i];
		for(mu=0; mu<2; mu++) overheatbath_single_link_U(conf, geo, param, n, mu, rng_state);
		overheatbath_single_site_z(conf, geo, param, n, rng_state);
	}
}	

void overheatbath_sweep_lattice(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, RNG_Param *rng_state)
{
	long i;
	int mu;
	for(i=0; i<param->d_volume ;i++)
	{
		for(mu=0; mu<2; mu++) overheatbath_single_link_U(conf, geo, param, i, mu, rng_state);
		overheatbath_single_site_z(conf, geo, param, i, rng_state);
	}
}

// over-heat-bath update of U and z fields
// Taken from "Monte Carlo simulation of CP^{N-1} models", M. Campostrini, P. Rossi, E. Vicari, Phys. Rev. D 46 (1992) 2647-2662

// over-heat-bath update of field U on link (i,mu)
void overheatbath_single_link_U(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, long const i, int const mu, RNG_Param *rng_state) 
{ 
	double a, b, b1, theta_new, cos_theta_new, sin_theta_new, sin_2_theta_old, cos_theta_old, mod_F_U;
	cmplx F_U;

	F_U=force_U(conf, geo, param, i, mu);
	mod_F_U = cmplx_abs(F_U); // |F_U|
	a = 2.0 * param->d_beta * ((double)N) * mod_F_U;

	cos_theta_old = creal( conj(conf->U[i][mu]) * F_U )/mod_F_U; // project U along F_U on link (i,mu)

	if (cos_theta_old > (1.0 - 1000.0*NUM_EPS) || mod_F_U==0) cos_theta_old = 1.0 - (1000.0*NUM_EPS); // avoid cos(theta_old) too close to 1 
	sin_2_theta_old=1.0-cos_theta_old*cos_theta_old; // this is sin^2(theta_old)
	if( sin_2_theta_old < NUM_EPS || mod_F_U==0) sin_2_theta_old = NUM_EPS; // avoid sin(theta_old) too close to 0

	// extract new polar angle between U and F_U
	theta_new=rng_theta_U(rng_state, a);      
	cos_theta_new=cos(theta_new); 
	sin_theta_new=sin(theta_new); 

	b  = sin_theta_new/sqrt(sin_2_theta_old); // sin(theta_new)/sin(theta_old)
	b1 = (cos_theta_new + cos_theta_old * b) / mod_F_U;

	conf->U[i][mu] *= -b; // U_new = -b U_old if F_U=0
	if(mod_F_U!=0) conf->U[i][mu] += b1*F_U; // U_new = -b U_old + b1 F_U if F_U!=0
	conf->U[i][mu] /= cmplx_abs(conf->U[i][mu]); // normalize U_new
}

// over-heat-bath update of field z on site i
void overheatbath_single_site_z(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, long const i, RNG_Param *rng_state) 
{
	double a, theta_new, cos_theta_new, sin_2_theta_new, sin_2_theta_old, cos_theta_old, b, b1, mod_F_z; 
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx a0;

	force_z(conf, geo, i, F_z);
	mod_F_z = vector_abs(F_z); // |F_z|
	a0=vector_scalar_product(conf->z[i], F_z); // project z along F_z on site i
	cos_theta_old = creal(a0) / mod_F_z; // Re( \conj(z) F_z )/|F_z| = cos(theta_old)

	sin_2_theta_old = 1.0-cos_theta_old*cos_theta_old; // this is sin^2(theta_old)
	a = 2.0 * param->d_beta * ((double)N) * mod_F_z; 

	theta_new=rng_theta_z(rng_state, a);
	cos_theta_new = cos(theta_new); 
	sin_2_theta_new = 1.0-cos_theta_new*cos_theta_new; // this is sin^2(theta_new)

	if ( sin_2_theta_old < NUM_EPS || mod_F_z==0) sin_2_theta_old = NUM_EPS; // avoid sin(theta_old) too close to 0
	if ( sin_2_theta_new < NUM_EPS ) sin_2_theta_new = NUM_EPS; // avoid sin(theta_new) too close to 0
	b = sqrt( sin_2_theta_new / sin_2_theta_old ); // sin(theta_new)/sin(theta_old)
	b1 = (cos_theta_new + cos_theta_old * b)/mod_F_z ;

	if (mod_F_z!=0) vector_linear_combination_real_coeff(conf->z[i], F_z, -b, b1); // z_new = -b z_old + b1 F_z if F_z!=0
	else vector_times_real_const(conf->z[i], -b); // z_new = -b z_old if F_z=0
	vector_normalization(conf->z[i]); // normalize z_new
}

// Metropolis test with Metropolis probability p. Returned values: 1=accepted, 0=rejected
int metropolis_test(RNG_Param *rng_state, double const p)
{
	double random_number;
	if (p>1) return 1;
	else
	{
		random_number=rand_num(rng_state);
		if (random_number<p) return 1;
		else return 0;
	}
}

// perform swaps of all couples od adjacent replicas
void swap(CPN_Conf * conf, CPN_Param const * const param, Acc_Swap *swap_counter, RNG_Param *rng_state)
{
	int label,r,s;
	for(label=0;label<((param->d_N_replica_pt)-1);label++)
	{
		r=label;
		s=label+1;
		swap_conf(conf, param, r, s, swap_counter, rng_state); // swp configurations with labels r and s
	}
}

// swap replicas with labels r and s
void swap_conf(CPN_Conf * conf,  CPN_Param const * const param, int const r, int const s, Acc_Swap *swap_counter, RNG_Param *rng_state)
{
	long n,cart_coord[2],pos;
	double E_r=0, E_s=0, log_p, p=0.0;
	int acc;

	// compute the action difference between the two replicas just along the defect or in the near-neighbors of it  
	for(n=0;n<param->d_L_defect;n++)
	{
		E_r+=swap_energy_contribution(&(conf[r]), param, n);
		E_s+=swap_energy_contribution(&(conf[s]), param, n);
	}

	// compute the boundary condition parameter conf[r].C[pos][0] in position pos=(L_0-1, 0) (could be any other position on the defect)
	// both contributions to the action differences (coming from the two terms of the action) have this constant in front
	cart_coord[0]=(param->d_size[0]-1);
	cart_coord[1]=0;
	pos=cart_to_si(cart_coord,param);

	log_p = (-2 * (param->d_beta) * ((double)N)) * (conf[r].C[pos][0]-conf[s].C[pos][0]) * (E_r - E_s); // log of the metropolis swap probability p
	p=exp(log_p); // metropolis swap probability p

	// Metropolis step for the swap of configurations
	swap_counter->num_swap[r]++; // increase the counter of the number of swaps tried during the simulation for the configuration with label r
	
	// Metropolis test
	acc=metropolis_test(rng_state, p);

	if(acc==1) // do the swap
	{
		// swap of configurations
		cmplx **z_aux, **U_aux;
		z_aux=conf[r].z;
		U_aux=conf[r].U;
		conf[r].z=conf[s].z;
		conf[r].U=conf[s].U;
		conf[s].z=z_aux;
		conf[s].U=U_aux;
		swap_counter->num_accepted_swap[r]++; // increase the counter of the number of swaps of the configuration with label r

		// swap of labels
		int aux_label;
		aux_label=conf[r].conf_label;
		conf[r].conf_label=conf[s].conf_label;
		conf[s].conf_label=aux_label;
	}
}
	
// compute contributions to energy difference for swap probability and site x = (x_0=L-1, x_1=n) (L = d_size[0]) and for links along mu=0
// generalized from "Fighting topological freezing in the two-dimensional CP^(N-1) model", M. Hasenusch, Phys. Rev. D 96 (2017) 5, 054504, 1706.04443 [hep-lat]
double swap_energy_contribution(CPN_Conf const * const conf, CPN_Param const * const param, long const n)
{
	cmplx a, b, c;
	double res;
	long pos1,pos2,pos3,cartcoord[2];

	// contribution of the first term in the Symanzik action

	// lexicographic coordinate of the site (x_0=0, x_1=n) 
	cartcoord[0]=0;
	cartcoord[1]=n;
	pos1=cart_to_si(cartcoord,param);

	// lexicographic coordinate of the site (x_0=L-1, x_1=n)
	cartcoord[0]=param->d_size[0]-1;
	cartcoord[1]=n;
	pos2=cart_to_si(cartcoord,param);

	a=vector_scalar_product(conf->z[pos1], conf->z[pos2]);
	a*=conj(conf->U[pos2][0]);

	// contribution of the second term in the Symanzik action

	// lexicographic coordinate of the site (x_0=1, x_1=n)
	cartcoord[0]=1;
	cartcoord[1]=n;
	pos3=cart_to_si(cartcoord,param);			

	// first addend			
	b=vector_scalar_product(conf->z[pos3], conf->z[pos2]);
	b*=conj(conf->U[pos1][0]) * conj(conf->U[pos2][0]);

	// lexicographic coordinate of the site (x_0=L-2, x_1=n)
	cartcoord[0]=param->d_size[0]-2;
	cartcoord[1]=n;
	pos3=cart_to_si(cartcoord,param);			

	// second addend		
	c=vector_scalar_product(conf->z[ pos1], conf->z[pos3]);
	c*=conj(conf->U[pos2][0]) * conj(conf->U[pos3][0]);

	// total
	res = c1 * creal(a) + c2 * creal(b) + c2 * creal(c);
	return res;			
}
	
// translation of one lattice spacing of the periodic configuration
// direction is chosen randomly, verse is always positive
void conf_translation(CPN_Conf *conf, Geometry const * const geo, CPN_Param const * const param, CPN_Conf *aux_conf, RNG_Param *rng_state)
{
	double aux;
	int dir=0,i,mu;
	long r;

	// extract random direction
	aux=rand_num(rng_state);
	aux=2.0*aux;
	for(i=0;i<2;i++)
	{
		if ( (aux>((double)i)) && (aux<((double)(i+1))) ) dir=i;
	}

	// aux_conf = conf
	copyconf(conf, param, aux_conf);

	// translation in direction +dir
	for(r=0;r<param->d_volume;r++)
	{
		vector_equal(conf->z[r], aux_conf->z[geo->up[r][dir]]);
		for(mu=0;mu<2;mu++) conf->U[r][mu]=aux_conf->U[ geo->up[r][dir] ][mu];
	}
}
  
void print_replicas_labels(FILE *fp, CPN_Conf const * const conf, CPN_Param const * const param)
{
	int r;
	if ( param->d_N_replica_pt == 1 )
	{
		(void) fp;   // to suppress compiler warnings about unused variables
		(void) conf; // to suppress compiler warnings about unused variables
	}
	else
	{
		fprintf(fp, "%ld      ",conf[0].update_index);
		for(r=0;r<(param->d_N_replica_pt);r++)
		{	  
			fprintf(fp,"%d ",conf[r].conf_label);
		}
	fprintf(fp,"\n");
	}
}

void init_swap_acceptances(Acc_Swap *swap_counter, CPN_Param const * const param)
{
	if(param->d_N_replica_pt==1)
	{
		swap_counter->num_accepted_swap=NULL;
		swap_counter->num_swap=NULL;
	}
	else
	{
		int i,err;
		err=posix_memalign( (void **) &(swap_counter->num_accepted_swap), (size_t) INT_ALIGN, (size_t) ((param->d_N_replica_pt)-1) * sizeof(long));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the acceptances array\n");
			exit(EXIT_FAILURE);
		}

		err=posix_memalign( (void **) &(swap_counter->num_swap), (size_t) INT_ALIGN, (size_t) ((param->d_N_replica_pt)-1) * sizeof(long));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the acceptances array\n");
			exit(EXIT_FAILURE);
		}

		for(i=0;i<((param->d_N_replica_pt)-1);i++) 
		{
			swap_counter->num_accepted_swap[i]=0;
			swap_counter->num_swap[i]=0;
		}
	}
}
  
void print_swap_acceptances(Acc_Swap const * const swap_counter, CPN_Param const * const param)
{
	if(param->d_N_replica_pt==1)
	{
		(void) swap_counter; // to suppress compiler warning of unused variable
	}
	else
	{
		FILE *fp;	
		double acc;
		int r;

		fp=fopen(param->d_swap_accept_file, "w");
		fprintf(fp, "#swap_from    swap_to    acceptance(%%)    swap_accepted    swap_tried\n");
		for(r=0;r<((param->d_N_replica_pt)-1);r++)
		{
			if(swap_counter->num_swap[r]!=0)
				acc=((((double)(swap_counter->num_accepted_swap[r]))/((double)(swap_counter->num_swap[r])))*100.0);
			else
				acc=0.0;
			fprintf(fp,"%d    %d    %f    %ld    %ld\n",r,r+1,acc,swap_counter->num_accepted_swap[r],swap_counter->num_swap[r]);
		}
		fclose(fp);	  
	}
}
	
void free_swap_acceptances(Acc_Swap *swap_counter, CPN_Param const * const param)
{
	if(param->d_N_replica_pt>1)
	{
		free(swap_counter->num_accepted_swap);
		free(swap_counter->num_swap);
	}
	else
	{
		(void) swap_counter; // to suppress compiler warning of unused variable
	}
}

void init_swap_track_file(FILE **swaptrackfilep, CPN_Param const * const param)
{    
	if (param->d_N_replica_pt==1)
	{
		(void) swaptrackfilep; // to suppress compiler warnings of unused variables
	}
	else
	{
		*swaptrackfilep=fopen(param->d_swap_tracking_file, "r");
		if(*swaptrackfilep!=NULL) // file exists
		{
			fclose(*swaptrackfilep);
			*swaptrackfilep=fopen(param->d_swap_tracking_file, "a");
		}
		else // file doesn't exist
		{
			*swaptrackfilep=fopen(param->d_swap_tracking_file, "w");
			fprintf(*swaptrackfilep, "# MC_step    conf_labels\n");
			fflush(*swaptrackfilep);
		}
		fflush(*swaptrackfilep);
	}
}

void end_swap_track_file(FILE **swaptrackfilep, CPN_Param const * const param)
{
	if ( param->d_N_replica_pt == 1 )
	{
		(void) swaptrackfilep; // to avoid compiler warning about unused variable
	}
	else
	{
		fclose(*swaptrackfilep);
	}
}

#endif
