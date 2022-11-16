#ifndef CPN_PARAM_H
#define CPN_PARAM_H

#include "macro.h"
#include "endianness.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// cpn param struct, stores all theory and simulation parameters
typedef struct CPN_Param {

	// lattice size
	int d_size[2];	// size[0] = time direction, size[1] = space direction
	
	// theory parameters
	double d_beta;
	double d_theta;
	
	// parallel tempering parameters
	int d_L_defect;			// defect length along the x direction;
	int d_N_replica_pt;		// numbers of replica used in parallel tempering
	
	// hierarchical update (parallel tempering)
	int d_N_hierarc_levels;	// number of hierarchical levels
	int *d_L_rect;			// d_L_rect is a vector of length d_N_hierarc_levels
							// d_L_rect[i] is the extension of the rectangle at the i-th hierarchical level
	int *d_N_sweep_rect;	// d_N_sweep_rect is vector of length d_N_hierarch_levels
							// d_N_sweep_rect[i] is the number of updating sweeps of the rectangle for the i-th hierarchical level
							
	// simulation parameters
	int d_MC_step;			// number of updating steps
	int d_measevery;		// number of updating steps between two measures
	int d_num_micro;		// number of over-relaxionation updating steps (microcanincal steps) for every over-heat-bath updating step
	int d_num_norm;			// normalize the configuration every d_num_norm updating steps
	
	// initialization and saving
	int d_start;					// initialize conf: 0 random cold, 1 random hot, 2 read conf from file
	int d_saveconf_backup_every;	// save current conf for backup every d_saveconf_backup_every updating steps
	
	// cooling
	int d_coolsteps;	// how many cooling steps are performed in total
	int d_coolevery;	// measure topological observables every <d_coolevery> cooling steps
	
	// random number generator parameters
	long d_seed;						// store rng seed
	int  d_rng_start;					// initialize rng state: 0 std initalization, 1 read rng state from file
	char d_rng_file[STD_STRING_LENGTH];	// rng state file name
	
	// topo potential parameters (multicanonic)
	double d_grid_step;									// topo potential is defined on a grid with spacing <grid_step>
	double d_grid_max;									// grid goes from -<grid_max> to <grid_max>
	char   d_topo_potential_file[STD_STRING_LENGTH];	// topo potential file
	int    d_stoc_single_site_upd;						// number of stochastic single site updates

	
	// output file names
	char d_conf_file[STD_STRING_LENGTH];				// conf file name
	char d_data_file[STD_STRING_LENGTH];				// non-topo data file name
	char d_topo_file[STD_STRING_LENGTH];				// topo data file name
	char d_log_file[STD_STRING_LENGTH];					// log file name
	char d_swap_accept_file[STD_STRING_LENGTH];			// swap acceptances file
	char d_swap_tracking_file[STD_STRING_LENGTH];		// swap history file
	char d_multicanonic_acc_file[STD_STRING_LENGTH];	// multicanonical acceptances file	
	
	// derived constants
	long d_volume;	// total lattice volume
	int  d_n_grid;	// total grid points

} CPN_Param;

// in lib/cpn_param.c
void remove_white_lines_and_comments(FILE *);
void read_input(char const * const, CPN_Param *);
void init_derived_constants(CPN_Param *);
void free_param(CPN_Param *);
void init_data_file(FILE **, CPN_Param const * const);
void init_topo_file(FILE **, CPN_Param const * const);
void print_simulation_details_cpn(char const * const, CPN_Param const * const, time_t const * const, time_t const * const, clock_t const, clock_t const);
void print_simulation_details_multicanonic_cpn(char const * const, CPN_Param const * const, time_t const * const, time_t const * const, clock_t const, clock_t const);

#endif
