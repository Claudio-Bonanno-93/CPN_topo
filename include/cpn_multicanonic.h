#ifndef CPN_MULTICANONIC_H
#define CPN_MULTICANONIC_H

#include "cpn_conf.h"

// in lib/cpn_update_multicanonic.c
void read_topo_potential(double **, CPN_Param const * const );
double compute_topo_potential(double const, double const * const, CPN_Param const * const);
void init_multicanonic_acc_file(FILE **, CPN_Param const * const );
void print_multicanonic_acceptance(CPN_Conf const * const, CPN_Param const * const, int const, long const, long const, double const, FILE *);
void multicanonic_parallel_tempering_with_hierarchic_update(CPN_Conf *, Rectangle const * const, Acc_Swap *, CPN_Param const * const, Geometry const * const,
                                                            CPN_Conf *, double const * const, FILE *, RNG_Param *);
void update_with_defect_stochastic(CPN_Conf * , Geometry const * const, CPN_Param const * const, double const * const, long *, long *, RNG_Param *);
void hierarchic_update_rectangle_stochastic(CPN_Conf *, Geometry const * const, CPN_Param const * const, int const,	int const,
                                       Rectangle const * const, Acc_Swap *, CPN_Conf *, double const * const, long *, long *, RNG_Param *);
void update_rectangle_stochastic(CPN_Conf *, Geometry const * const, CPN_Param const * const, Rectangle const * const,
                                 double const * const, long *, long *, RNG_Param *);
void init_topo_charge(CPN_Conf * , Geometry const * const, CPN_Param const * const);
void refresh_topo_charge_replica(CPN_Conf *, Geometry const * const, CPN_Param const * const);
void refresh_topo_charge(CPN_Conf *, Geometry const * const, CPN_Param const * const);
double delta_Q(CPN_Conf const * const, Geometry const * const, cmplx const, long const, int const);
double metropolis_prob_multicanonic(double const, double const, CPN_Param const * const, double const * const);
int multicanonic_update_single_link_U(CPN_Conf *, Geometry const * const, CPN_Param const * const, double const * const,
										long const, int const, RNG_Param *, double const);

#endif
