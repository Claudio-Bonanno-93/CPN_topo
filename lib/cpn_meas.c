#ifndef CPN_MEAS_C
#define CPN_MEAS_C

#include "../include/cpn_conf.h"
#include <time.h>

// measure local observables on the periodic copy
void perform_measures_localobs(CPN_Conf * conf, Geometry const * const geo, 
                               CPN_Param const * const param, FILE *datafilep, FILE *topofilep, CPN_Conf *aux_conf)
{
	int cool_step=0, i;
	double magn_susc[2], Q[3], chi_p[3];
	double energy;

	energy=energy_density(conf, geo, param);
	magnetic_susceptibility(conf, param, magn_susc);

	fprintf(datafilep, "%ld %.16lf %.16lf %.16lf\n", conf->update_index, energy, magn_susc[0], magn_susc[1]);
	fflush(datafilep);

	// compute topological observables of hot configuration
	for (i=0; i<3; i++) // 1=geometric charge U, 2=geometric charge z, 3=non-geometric plaquette charge
	{
	Q[i] = topo_charge(conf, geo, param, i); // compute topological charge using i^th discretization
	chi_p[i] = chi_prime(conf, geo, param, i); // compute chi' using i^th discretization
	}	
	// print topological observable of hot configuration
	fprintf(topofilep, "%ld %d", conf->update_index, cool_step);
	for (i=0; i<3; i++) fprintf(topofilep, " %.16lf", Q[i]);
	for (i=0; i<3; i++) fprintf(topofilep, " %.16lf", chi_p[i]);
	fprintf(topofilep, "\n");

	// refresh stored topo charge of periodic configuration (used only for multicanonic)
	conf->stored_topo_charge = Q[0];

	// aux_conf = conf (used for cooling)
	copyconf(conf, param, aux_conf);

	// perform cooling on aux_conf
	for (cool_step=1; cool_step<(param->d_coolsteps+1); cool_step++) 
    { 
		cooling(aux_conf, geo, param); // perform 1 cooling step
		if ( ( cool_step % param->d_coolevery ) == 0 ) // perform measures on cooled conf
		{ 
			// compute topological observables of cooled configuration
			for (i=0; i<3; i++)
			{
				Q[i] = topo_charge(aux_conf, geo, param, i);
				chi_p[i] = chi_prime(aux_conf, geo, param, i);
			}
			// print topological observable of cooled configuration
			fprintf(topofilep, "%ld %d", conf->update_index, cool_step);
			for (i=0; i<3; i++) fprintf(topofilep, " %.16lf", Q[i]);
			for (i=0; i<3; i++) fprintf(topofilep, " %.16lf", chi_p[i]);
			fprintf(topofilep, "\n");
		}
	}

	fflush(topofilep);
}

// compute plaquette Pi_{mu nu}(i) on site i and plane (mu,nu)
// Pi_{mu nu}(i) = U(i)_mu U(i+mu)_nu conj( U(i+nu) )_mu conj( U(i) )_nu
cmplx plaquette(CPN_Conf const * const conf, Geometry const * const geo, long const i, int const mu, int const nu)
{
	return( conf->U[i][mu] * conf->U[geo->up[i][mu]][nu] * conj(conf->U[geo->up[i][nu]][mu]) * conj(conf->U[i][nu]) );
}

// compute the energy density for single link E = S_{Symanzik}(theta=0) / (2 V N beta)
double energy_density(CPN_Conf const * const conf, Geometry const * const geo, CPN_Param const * const param) 
{ 
	int i, mu;
	double e = 0.0;
	cmplx e1 = 0.0 + I * 0.0, e2 = 0.0 + I * 0.0; 
	cmplx aux1, aux2;

	for (i=0; i<param->d_volume; i++) 
	{ 
		for(mu=0 ; mu<2 ; mu++)
		{
			aux1=vector_scalar_product(conf->z[geo->up[i][mu]], conf->z[i]); // conj( z(i+mu) ) * z(i)
			e1 += conj(conf->U[i][mu]) * aux1;

			aux2=vector_scalar_product(conf->z[geo->up[geo->up[i][mu]][mu]], conf->z[i]); // conj( z(i+2mu) ) * z(i)
			e2 += conj(conf->U[geo->up[i][mu]][mu] * conf->U[i][mu]) * aux2; 
		}
	}

	e = 2.0*(c1*creal(e1) + c2*creal(e2)); 
	e /= (double)(param->d_volume);
	e = 2.0*(c1 + c2) - e/2.0;
	return e; 
}

// compute the geometric topological charge density expressed in terms of the scalar field z on site i
double geo_topo_charge_z_density(CPN_Conf const * const conf, Geometry const * const geo, long const i)
{
	cmplx aux_1, aux_2, aux_3, p1, p2;
	double q; 
	int mu=0;
	int nu=1-mu; 

	aux_1 = vector_scalar_product(conf->z[geo->up[ geo->up[i][mu]][nu]], conf->z[i]);
	aux_2 = vector_scalar_product(conf->z[geo->up[i][mu]], conf->z[geo->up[geo->up[i][mu]][nu] ]);
	aux_3 = vector_scalar_product(conf->z[i], conf->z[geo->up[i][mu]]);
	p1 = aux_1 * aux_2 * aux_3;

	aux_1 = vector_scalar_product(conf->z[geo->up[i][nu]], conf->z[i]);
	aux_2 = vector_scalar_product(conf->z[geo->up[geo->up[i][mu]][nu]], conf->z[geo->up[i][nu]]);
	aux_3 = vector_scalar_product(conf->z[i], conf->z[geo->up[geo->up[i][mu]][nu]]);
	p2 = aux_1 * aux_2 * aux_3;

	q = -1.0 * ( arg(p1) + arg(p2) ) / (2.0*pi);
	return q;
}

// compute the geometric topological charge expressed in terms of the gauge field U on site i
double geo_topo_charge_U_density(CPN_Conf const * const conf, Geometry const * const geo, long const i)
{
	int mu=0;
	int nu=1-mu;
	cmplx plaq;
	double q;

	plaq = plaquette(conf, geo, i, mu, nu); // Plaq = plaquette(i)_{01}
	q = arg(plaq) / (2.0*pi); // {Im log(plaq) } / (2 pi)
	return q;
}

// compute the non-geometric topological charge expressed in terms of the plaquettes on site i
double plaq_topo_charge_density(CPN_Conf const * const conf, Geometry const * const geo, long const i)
{
	int mu=0;
	int nu=1-mu;
	double q;

	q = cimag( plaquette(conf, geo, i, mu, nu) ); // Im{ plaquette(i)_{01} } 
	q /= (2.0*pi); 
	return q;
}

// compute the total topological charge
double topo_charge(CPN_Conf const * const conf, Geometry const * const geo, CPN_Param const * const param, int const which_charge)
{
	// pointer to desired discretization of the topological charge density
	double (*q_ptr)(CPN_Conf const * const, Geometry const * const, long const);
	if (which_charge==0) {q_ptr = &geo_topo_charge_U_density;}
	if (which_charge==1) {q_ptr = &geo_topo_charge_z_density;}
	if (which_charge==2) {q_ptr = &plaq_topo_charge_density;}

	double Q = 0.0; 
	long i; 
	for(i=0; i<param->d_volume; i++) {Q += ((*q_ptr)(conf, geo, i));} // Q = sum_i { q(i) }
	return Q;
}

// compute quantity G = (1/4) sum_{x} |x|^2 q(x)q(0) needed to compute chi'
// the mean value of this sum in the continuum limit is chi' = (1/4) \int d^2x |x|^2 <q(x)q(0)> and is related to the first derivative with respect to k=p^2 of
// the Fourier transform FG of the topological charge density correlator: chi' = - (1/4) lim_{k->0} d FG(k) / dk
// where FG(p^2) = \int d^2 x e^(ipx) <q(x)q(0)>
double chi_prime(CPN_Conf const * const conf, Geometry const * const geo, CPN_Param const * const param, int const which_charge)
{
	// pointer to desired discretization of the topological charge density
	double (*q_ptr)(CPN_Conf const * const, Geometry const * const, long const);
	if (which_charge==0) {q_ptr = &geo_topo_charge_U_density;}
	if (which_charge==1) {q_ptr = &geo_topo_charge_z_density;}
	if (which_charge==2) {q_ptr = &plaq_topo_charge_density;}

	double d2, G = 0.0;
	long i;
	for (i=0; i<param->d_volume; i++)
	{
		d2 = square_distance(i, 0, param); // i->x ==> d2 = d(x,0)^2 = |x|^2 = |i|^2 where distance is computed on the torus
		G += ((*q_ptr)(conf, geo, i))*d2; // sum_i { q(i) |i|^2 }
	}
	G *= ((*q_ptr)(conf, geo, 0));  // sum_i { q(i)q(0) |i|^2 }
	G /= 4.0; // G = (1/4) sum_i { q(i)q(0) |i|^2 }
	return G;
}

// compute quantities needed for the "magnetic susceptibility" at momentum p=0 (stored in magn_susc[0]) and at momentum p = q = 2pi/L (stored in magn_susc[1])
// q is the smallest momentum possible on a lattice with L sites, where L is the smallest size of the lattice
// the magnetic susceptibility can be expressed in terms of G(p) = Fourier transform of G(x) = P_(ij)(x) P_(ij)(0) - 1/N, where P_(ij)(x) = z*_i(x) z_j(x)
// magnetic susceptibility at p=0: chi_m(0) = <G(p=0)> (this quantity is trivially real)
// magnetic susceptibility at p=q: chi_m(q) = <G(p=q)> (here I take just the real part, as the imaginary part averages to zero over the ensamble)
// mang_susc[0] = G(p=0), magn_susc[1]=G(p=q), the correlation length of the system can be expressed in terms of chi_m(q)/chi_m(0)
void magnetic_susceptibility(CPN_Conf const * const conf, CPN_Param const * const param, double * magn_susc)
{ 
	int i;
	double L, sum1=0.0, q;
	cmplx a;
	cmplx sum2 = 0.0 + I * 0.0;
	long x[2];

	// L = min(Lx,Lt)
	if (param->d_size[0] > param->d_size[1]) L = ( (double) param->d_size[1] ); 
	else L = ( (double) param->d_size[0] );
	q = 2.0*pi/L;

	for(i=0; i<param->d_volume; i++)
	{
		a=vector_scalar_product(conf->z[i], conf->z[0]);
		sum1 += cmplx_norm(a);		
		si_to_cart(x, i, param); // i->x
		sum2 += ( cmplx_norm(a) - 1.0/((double)N) )*cexp( I * ( q * ( (double)x[0] ) ));
	}
	sum1 -= ((double)(param->d_volume))/((double)N);

	magn_susc[0] = sum1;
	magn_susc[1] = creal(sum2);
}

// perform a single cooling step on the given conf (this function should receive an aux conf)
// cooling = each field is aligned to its local force: z = F_z / |F_z|, U = F_U / |F_U|
// forces are determined using the non-improved lattice action S_L \propto sum_{i, mu} U(i)_mu conj(z)(i+mu) z(i)
void cooling (CPN_Conf *conf, Geometry const * const geo, CPN_Param const * const param) 
{ 
	int i, mu;
	cmplx F_U; 
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx aux[N] __attribute__((aligned(DOUBLE_ALIGN))); 

	// cool U fields
	for (i=0; i<param->d_volume ; i++) 
	{
		for(mu=0 ; mu<2 ; mu++)
		{
		F_U=vector_scalar_product(conf->z[geo->up[i][mu]], conf->z[i]); // F_U = conj(z(i+mu)) z(i)
		conf->U[i][mu] = F_U / cmplx_abs(F_U); // align U along the local force on link (i,mu): U = F_U/|F_U|
		}
	}

	// cool z fields
	for (i=0; i<param->d_volume; i++)
	{
		vector_zero(F_z); // F_z = 0
		for (mu=0; mu<2; mu++)
		{
			// aux = U(i)_mu * z(i+mu) + conj( U(i-mu)_mu ) * z(i-mu)
			vector_linear_combination_cmplx_coeff(aux, conf->z[geo->up[i][mu]], conf->z[geo->dn[i][mu]], conf->U[i][mu], conj(conf->U[geo->dn[i][mu]][mu]));
			vector_sum(F_z, aux); // F_z += aux;
		}
		vector_normalization(F_z); // F_z -> F_z/|F_z|
		vector_equal(conf->z[i], F_z); // align z along the local force on site i: z = F_z/|F_z|
	}
}

#endif
