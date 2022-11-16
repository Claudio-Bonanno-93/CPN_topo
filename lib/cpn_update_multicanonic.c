#ifndef CPN_UPDATE_MULTICANONIC_C
#define CPN_UPDATE_MULTICANONIC_C

#include "../include/cpn_multicanonic.h"

// read topo potential from file
void read_topo_potential(double **grid, CPN_Param const * const param)
{
	int i, j, err;
	double x, V;
	FILE *fp;

	fp=fopen(param->d_topo_potential_file, "r");
	if( fp==NULL )
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_topo_potential_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}	
	
	err=posix_memalign((void **) grid, (size_t) DOUBLE_ALIGN, (size_t) param->d_n_grid * sizeof(double));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the grid array!\n");
		exit(EXIT_FAILURE);
	}

	// read x and V(x) from topo_potential file
	for (i=0; i<param->d_n_grid; i++)
	{
		err=fscanf(fp, "%lf %lf", &x, &V); 
		if(err!=2)
		{
			printf("Error: can't read the %d-th element of the file (%s, %d)\n", i,__FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		j=(int)floor((x+param->d_grid_max+(param->d_grid_step/2.0))/param->d_grid_step);
		if (i!= j)
		{
			printf("Error: found %d (%lf) when expecting %d (%s, %d)\n", j, x, i, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		(*grid)[i]=V;
	}
	fclose(fp);
}

// compute topo potential in point x
double compute_topo_potential(double const x, double const * const grid, CPN_Param const * const param)
{	  
	int i_grid;
	double x0, m, q;

	// find index of nearest grid point to x
	i_grid=(int)(floor((x+param->d_grid_max)/param->d_grid_step));

	if(i_grid>=0 && i_grid<param->d_n_grid) // if x inside the barriers compute V(x) with a linear interpolation
	{
		// perform linear interpolation: V(x) = V(x0) + [dV/dx|(x0)] (x-x0)
		x0=i_grid*param->d_grid_step-param->d_grid_max;
		m=(grid[i_grid+1]-grid[i_grid])/param->d_grid_step; // dV/dx|(x0) = [ V(x0+step) - V(x0) ] / step
		q=grid[i_grid]-m*x0; // V(x0) - [dV/dx|(x0)] x0
		return q+m*x; // V(x0) + [dV/dx|(x0)] (x-x0) 
	}
	else // if x outside the barriers just saturate to extreme values
	{
		if(i_grid<0) return grid[0];
		else return grid[param->d_n_grid-1];
	}
}

// initialize multicanonic acceptances file
void init_multicanonic_acc_file(FILE **filep, CPN_Param const * const param)
{
	int i;
	
	(*filep)=fopen(param->d_multicanonic_acc_file, "r");
	if((*filep)!=NULL) // file exists
	{
		fclose((*filep));
		(*filep)=fopen(param->d_multicanonic_acc_file, "a");
	}
	else // file doesn't exist, write first line
	{
		(*filep)=fopen(param->d_multicanonic_acc_file, "w");
		fprintf((*filep), "# %f ", param->d_beta);
		for(i=0; i<2; i++) fprintf(*filep, "%d ", param->d_size[i]);
		fprintf((*filep), "\n");
	}
	fflush(*filep);
}

// print metropolis acceptance of multicanonic algorithm on file
void print_multicanonic_acceptance(CPN_Conf const * const conf, CPN_Param const * const param, int const replica_index,
                                   long const num_accept, long const num_metropolis, double const mean_acc, FILE *multicanonic_acc_filep)
{
	if( conf->update_index % param-> d_measevery == 0)
	{
		if ( param->d_N_replica_pt == 1)
		{
			fprintf(multicanonic_acc_filep, "%ld %ld %ld %.12lf\n", conf->update_index, num_accept, num_metropolis, mean_acc);
			fflush(multicanonic_acc_filep);
		}
		else
		{
			fprintf(multicanonic_acc_filep, "%ld %d %ld %ld %.12lf\n", conf->update_index, replica_index, num_accept, num_metropolis, mean_acc);
			fflush(multicanonic_acc_filep);
		}
	}
}

// perform a single step of parallel tempering with hierarchic update with a multicanonical approach
void multicanonic_parallel_tempering_with_hierarchic_update(CPN_Conf *conf, Rectangle const * const most_update, Acc_Swap *swap_counter, CPN_Param const * const param,
                                                    Geometry const * const geo, CPN_Conf *aux_conf, double const * const grid, FILE *acc_filep, RNG_Param *rng_state)
{
	int i, begin_hierarc_level=0;
	long num_accept=0, num_metropolis=0;
	double mean_acc;

	// iterate over replicas
	for(i=0; i<param->d_N_replica_pt; i++)
	{
		/*----- Parallel Tempering updating step for i-th replica with multicanonic -----*/
		update_with_defect_stochastic(&(conf[i]), geo, param, grid, &num_accept, &num_metropolis, rng_state);	
		if(param->d_N_replica_pt>1)
		{
			swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
			if ( param->d_N_hierarc_levels > 0 )
			{
				hierarchic_update_rectangle_stochastic(conf, geo, param, i, begin_hierarc_level, most_update, swap_counter,
																								aux_conf,grid, &num_accept, &num_metropolis, rng_state);
			}
		}
		/*-------------------------------------------------------------------------------*/
		
		// print mean multicanonic acceptance over a single update step for i-th replica
		if (num_metropolis != 0) mean_acc=( (double) num_accept)/( (double) num_metropolis);
		else mean_acc=0.0;
		print_multicanonic_acceptance(conf, param, i, num_accept, num_metropolis, mean_acc, acc_filep);
		
		// increase counter of update steps for i-th replica
		conf[i].update_index++;
	}
}


// perform a stochastic update
// 1) random choice of update type: over-relaxation with probability p and over-heat-bath with probability 1-p, p=num_micro/(num_micro + 1)
// 2) random uniform choice of site i
// 3) random uniform choice of the field to update: z, U_0 or U_1
void update_with_defect_stochastic(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, double const * const grid,
                                   long *num_accept, long *num_metropolis, RNG_Param *rng_state)
{
	int acc,mu;
	double x, y, prob_micro=((double)param->d_num_micro)/((double)( 1 + param->d_num_micro));
	long i, n;

	for(i=0; i<param->d_stoc_single_site_upd; i++)
	{
		acc=0;
		n=(long int) floor(rand_num(rng_state) * ((double)param->d_volume)); // random site
		x=rand_num(rng_state); // random number for field choice
		y=rand_num(rng_state); // random number for update choice

		if (x<(1.0/3.0)) // update z field
		{
			if (y < prob_micro) microcanonic_single_site_z(conf, geo, n);
			else overheatbath_single_site_z(conf, geo, param, n, rng_state);
		}	
		else // update U field
		{
			if(x > (1.0/3.0) && x < (2.0/3.0)) mu=0; // update U_0 field
			else mu=1; // update U_1 field
			acc=multicanonic_update_single_link_U(conf, geo, param, grid, n, mu, rng_state, y-prob_micro);
			(*num_metropolis)++;
		}
	(*num_accept)+=acc;
	}
}

// perform a hierarchic update on rectangles for conf[conf_label]
void hierarchic_update_rectangle_stochastic(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param, int const conf_label, int const hierarc_level,
									Rectangle const * const most_update, Acc_Swap *swap_counter, CPN_Conf *aux_conf,
									double const * const grid, long *num_accept, long *num_metropolis, RNG_Param *rng_state)
{
	int j;
	if(hierarc_level==((param->d_N_hierarc_levels)-1))
	{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++) 
		{
			update_rectangle_stochastic(&(conf[conf_label]), geo, param, &(most_update[hierarc_level]), grid, num_accept, num_metropolis, rng_state);
			if(param->d_N_replica_pt>1) swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
		}
	} // end if
	else
	{
		for(j=0;j<param->d_N_sweep_rect[hierarc_level];j++)
		{
			update_rectangle_stochastic(&(conf[conf_label]), geo, param, &(most_update[hierarc_level]), grid, num_accept, num_metropolis, rng_state);   
			if(param->d_N_replica_pt>1) swap(conf, param, swap_counter, rng_state);
			conf_translation(&(conf[0]), geo, param, aux_conf, rng_state);
			hierarchic_update_rectangle_stochastic(conf, geo, param, conf_label, (hierarc_level+1), most_update, swap_counter,
													aux_conf, grid, num_accept, num_metropolis, rng_state);
		}
	} // end else
}

void update_rectangle_stochastic(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param,
                                 Rectangle const * const most_update, double const * const grid, long *num_accept, long *num_metropolis, RNG_Param *rng_state)
{
	int acc, mu;
	double x, y, prob_micro=((double)param->d_num_micro)/((double)( 1 + param->d_num_micro));
	long i, n, m, num_rect_upd;

	if (param->d_stoc_single_site_upd < most_update->d_vol_rect) num_rect_upd=param->d_stoc_single_site_upd;
	else num_rect_upd=most_update->d_vol_rect;

	for(i=0; i<num_rect_upd; i++)
	{
		acc=0;
		m=(long int) floor(rand_num(rng_state) * ((double)most_update->d_vol_rect));
		n=most_update->rect_sites[m]; // random rectangle site
		x=rand_num(rng_state); // random number for field choice
		y=rand_num(rng_state); // random number for update choice

		if (x<(1.0/3.0)) // update z field
		{
			if(y < prob_micro) microcanonic_single_site_z(conf, geo, n); 
			else overheatbath_single_site_z(conf, geo, param, n, rng_state);			
		}
		else // update U field
		{
			if(x > (1.0/3.0) && x < (2.0/3.0)) mu=0; // update U_0 field
			else mu=1; // update U_1 field
			acc=multicanonic_update_single_link_U(conf, geo, param, grid, n, mu, rng_state, y-prob_micro);
			(*num_metropolis)++;
		}
		*num_accept=(*num_accept)+acc;
	}
}

// initialize topo_charge for all replicas
void init_topo_charge(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param)
{
	int i=0;
	for (i=0; i<param->d_N_replica_pt; i++) refresh_topo_charge(&(conf[i]), geo, param);
}

// refresh topo_charge for non-periodic replica
void refresh_topo_charge_replica(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param)
{
	// starts from 1 because the topo charge of replica with i=0 is refreshed before measuring observables
	int i=1;
	for (i=1; i<param->d_N_replica_pt; i++) refresh_topo_charge(&(conf[i]), geo, param);
}

// refresh topo_charge of given replica
void refresh_topo_charge(CPN_Conf *conf, Geometry const * const geo, CPN_Param const * const param)
{
	double Q;
	Q=topo_charge(conf, geo, param, 0);
	conf->stored_topo_charge=Q;
}

// compute the variation of the geometric topological charge when the link (i,mu) is updated from U_old to U_new
double delta_Q_upd(CPN_Conf const * const conf, Geometry const * const geo, cmplx const U_old, long const i, int const mu)
{
	double delta_Q, q_a, q_b, q_c, q_d;
	cmplx a, b, c, d, stpl_up, stpl_dn, U_new;

	U_new=conf->U[i][mu];
	
	stpl_up = staple_up(conf, geo, i, mu);
	stpl_dn = staple_down(conf, geo, i, mu);
	
	if (mu==0)
	{
		a = U_new*stpl_up;
		q_a = arg(a);

		b = conj(U_new*stpl_dn);
		q_b = arg(b);

		c = U_old*stpl_up;
		q_c = arg(c);

		d = conj(U_old*stpl_dn);
		q_d = arg(d);
	}
	else
	{
		a = conj(U_new*stpl_up);
		q_a = arg(a);

		b = U_new*stpl_dn;
		q_b = arg(b);

		c = conj(U_old*stpl_up);
		q_c = arg(c);

		d = U_old*stpl_dn;
		q_d = arg(d);
	}

	delta_Q = q_a + q_b - q_c - q_d;
	delta_Q /= (2.0*pi);
	return delta_Q;
}

// compute the multicanonic Metropolis probability p=exp(delta V) where V=topo potential
double metropolis_prob_multicanonic(double const Q_new, double const Q_old, CPN_Param const * const param, double const * const grid)
{
	double V_old, V_new, delta_V, p;

	V_old = compute_topo_potential(Q_old, grid, param);
	V_new = compute_topo_potential(Q_new, grid, param);
	delta_V = V_new-V_old;
	p=exp(-delta_V);
	return p;
}

// stochastic over-relaxation / over-heat-bath update of field U on link (i, mu) followed by a Metropolis test with p_metro=exp(- delta topo_potential)
int multicanonic_update_single_link_U(CPN_Conf * conf, Geometry const * const geo, CPN_Param const * const param,
                                            double const * const grid, long const i, int const mu, RNG_Param *rng_state, double const which_update)
{
	int acc;
	double Q_old, Q_new, p;
	cmplx U_old;

	// store topological charge and U on link (i,mu) before update
	Q_old=conf->stored_topo_charge;
	U_old=conf->U[i][mu];

	// perform update
	// which_update = y - prob_micro   ==>   which_update < 0 --> over-relaxation, which_update > 0 --> over-heat-bath
	if (which_update < 0) microcanonic_single_link_U(conf, geo, param, i, mu);
	else overheatbath_single_link_U(conf, geo, param, i, mu, rng_state);

	// perform multicanonic Metropolis test
	Q_new = Q_old + delta_Q_upd(conf, geo, U_old, i, mu);
	p=metropolis_prob_multicanonic(Q_new, Q_old, param, grid);
	acc=metropolis_test(rng_state, p);

	if (acc == 0) conf->U[i][mu] = U_old;	 // if Metropolis is refused go back to original link
	else conf->stored_topo_charge = Q_new; // if Metropolis is accepted store the new topological charge of the conf

	return acc;
}

#endif
