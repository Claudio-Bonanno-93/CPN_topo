#ifndef CPN_CONF_H
#define CPN_CONF_H

#include "geometry.h"
#include "cpn_cmplx_op.h"

#ifdef ENABLE_MD5_HASH
#include <openssl/md5.h>
#endif

// CPN conf struct, stores field configurations
typedef struct CPN_Conf {

	long update_index; // counts number of updates performed
	
	cmplx ** z, ** U;  // CP^{N-1} fields z[space-time volume][N] and U[space-time volume][space-time dimension=2]
	
	//for multicanonic
	double stored_topo_charge; // stores topological charge 
	
	// for parallel tempering
	double **C;       // C[space-time volume][space-time-dimension=2], this factor changes the boundary condition for the links crossing the defect
	int conf_label;   // stores the label of the replica to keep track of the swaps

} CPN_Conf;

// Hasenbusch parallel tempering swap struct, stores the accepted swaps and the total swaps proposed during Hasenbusch parallel tempering
typedef struct Acc_Swap {

  long *num_accepted_swap;
  long *num_swap;

} Acc_Swap;

// in lib/cpn_conf_def.c
void init_CPN_replicas(CPN_Conf **, CPN_Param const * const, RNG_Param *);
void allocate_CPN_conf(CPN_Conf *, CPN_Param const * const);
void init_CPN_conf(CPN_Conf *, CPN_Param const * const, char const * const, RNG_Param *); 
void init_bound_cond(CPN_Conf *, int const, CPN_Param const * const);
int is_on_defect(long const, CPN_Param const * const);
void normalize_replicas(CPN_Conf *, CPN_Param const * const); 
void normalize_CPN_conf(CPN_Conf *, CPN_Param const * const);
void copyconf(CPN_Conf const * const, CPN_Param const * const, CPN_Conf *);
void write_replicas(CPN_Conf const * const , CPN_Param const * const);
void write_replicas_backup(CPN_Conf const * const, CPN_Param const * const);
void write_CPN_conf_on_file(CPN_Conf const * const, CPN_Param const * const, char const * const);
void read_CPN_conf_from_file(CPN_Conf *, CPN_Param const * const, char const * const);
void compute_MD5_hash_conf(char *, CPN_Conf const * const, CPN_Param const * const);
void free_CPN_replicas(CPN_Conf *, CPN_Param const * const);
void free_bound_cond(CPN_Conf *, CPN_Param const * const);
void free_CPN_conf(CPN_Conf *, CPN_Param const * const);

// in lib/cpn_meas.c
void perform_measures_localobs(CPN_Conf  * , Geometry const * const, CPN_Param const * const, FILE *, FILE *, CPN_Conf *);
cmplx plaquette(CPN_Conf const * const, Geometry const * const, long const, int const mu, int const nu);
double energy_density(CPN_Conf const * const, Geometry const * const, CPN_Param const * const);
double geo_topo_charge_z_density(CPN_Conf const * const, Geometry const * const, long const);
double geo_topo_charge_U_density(CPN_Conf const * const, Geometry const * const, long const); 
double plaq_topo_charge_density(CPN_Conf const * const, Geometry const * const, long const);
double topo_charge(CPN_Conf const * const, Geometry const * const, CPN_Param const * const, int const);
double chi_prime(CPN_Conf const * const, Geometry const * const, CPN_Param const * const, int const);
void magnetic_susceptibility(CPN_Conf const * const, CPN_Param const * const, double *); 
void cooling(CPN_Conf *, Geometry const * const, CPN_Param const * const);

// in lib/cpn_update.c
void parallel_tempering_with_hierarchic_update(CPN_Conf *, Rectangle const * const, Acc_Swap *,
                                               CPN_Param const * const, Geometry const * const, CPN_Conf *, RNG_Param *);
void hierarchic_update_rectangle(CPN_Conf *, Geometry const * const, CPN_Param const * const, int const, int const,
                                 Rectangle const * const, Acc_Swap *, CPN_Conf *, RNG_Param *);
void update_rectangle(CPN_Conf *, Geometry const * const, CPN_Param const * const, Rectangle const * const, RNG_Param *);
void update_with_defect(CPN_Conf *, Geometry const * const, CPN_Param const * const, RNG_Param *);
cmplx staple_up(CPN_Conf const * const, Geometry const * const, long const, int);
cmplx staple_down(CPN_Conf const * const, Geometry const * const, long const, int);
cmplx force_U(CPN_Conf const * const, Geometry const * const, CPN_Param const * const, long const, int const); 
void force_z(CPN_Conf const * const, Geometry const * const, long const, cmplx *); 
void microcanonic_sweep_rectangle(CPN_Conf *, Geometry const * const, CPN_Param const * const, Rectangle const * const);
void microcanonic_sweep_lattice(CPN_Conf *, Geometry const * const, CPN_Param const * const);
void microcanonic_single_link_U(CPN_Conf *, Geometry const * const, CPN_Param const * const, long const, int const);
void microcanonic_single_site_z(CPN_Conf *, Geometry const * const, long const);
void overheatbath_sweep_rectangle(CPN_Conf *, Geometry const * const, CPN_Param const * const, Rectangle const * const, RNG_Param *);
void overheatbath_sweep_lattice(CPN_Conf *, Geometry const * const, CPN_Param const * const, RNG_Param *); 
void overheatbath_single_link_U(CPN_Conf *, Geometry const * const, CPN_Param const * const, long const, int const, RNG_Param *);
void overheatbath_single_site_z(CPN_Conf *, Geometry const * const, CPN_Param const * const, long const, RNG_Param *);
int metropolis_test(RNG_Param *, double const);
void swap(CPN_Conf *, CPN_Param const * const, Acc_Swap *, RNG_Param *);
void swap_conf(CPN_Conf *, CPN_Param const * const, int const, int const, Acc_Swap *, RNG_Param *);
double swap_energy_contribution(CPN_Conf const * const, CPN_Param const * const, long const);
void conf_translation(CPN_Conf *, Geometry const * const, CPN_Param const * const, CPN_Conf *, RNG_Param *);
void print_replicas_labels(FILE *, CPN_Conf const * const, CPN_Param const * const);
void init_swap_acceptances(Acc_Swap *, CPN_Param const * const);
void print_swap_acceptances(Acc_Swap const * const, CPN_Param const * const);
void free_swap_acceptances(Acc_Swap *, CPN_Param const * const);
void init_swap_track_file(FILE **, CPN_Param const * const);
void end_swap_track_file(FILE **, CPN_Param const * const);

#endif
