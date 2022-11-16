#ifndef CPN_C
#define CPN_C

#include "../include/macro.h" // NOTE: must be included before stdlib.h to activate posix_memalign correctly
#include "../include/cpn_conf.h"
#include "../include/cpn_cmplx_op.h"
#include "../include/cpn_param.h"
#include "../include/geometry.h"
#include "../include/rng.h"
#include "../include/endianness.h"

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

void real_main(char *input_file_name) 
{
	CPN_Conf *conf;
	CPN_Conf aux_conf;
	CPN_Param param;	
	Geometry geo;
	Rectangle *most_update;
	Acc_Swap swap_counter;
	RNG_Param rng_state;
	time_t start_date, finish_date;
	clock_t start_time, finish_time;
	FILE *datafilep, *swaptrackfilep, *topofilep;
	int i;

	// read input file
	read_input(input_file_name, &param);

	// initialization of rng state
	init_rng_state(&rng_state, &param);

	// open data file
	init_data_file(&datafilep, &param);

	// open topo data file
	init_topo_file(&topofilep, &param);

	// open swap tracking file
	init_swap_track_file(&swaptrackfilep, &param);

	// initialize geometry
	init_geometry(&geo, &param);

	// initialize lattice replicas
	init_CPN_replicas(&conf, &param, &rng_state);

	// initialize aux conf (will be used for cooling and for periodic conf translations)
	allocate_CPN_conf(&aux_conf, &param);

	// initialize rectangles for hierarchic updates
	init_rectangles_hierarchic_upd(&most_update, &param);

	// initialize swap acceptances array
	init_swap_acceptances(&swap_counter, &param);

 	// Monte Carlo begins
	time(&start_date);
	start_time=clock();
	for (i=0; i<param.d_MC_step; i++)  
	{
		// perform a single step of parallel tempering with hierarchic updates and print state of replica swap
		parallel_tempering_with_hierarchic_update(conf, most_update, &swap_counter, &param, &geo, &aux_conf, &rng_state);

		// normalize the lattice fields of all replicas
		if ( (conf[0].update_index) % param.d_num_norm == 0) normalize_replicas(conf,&param); 
		
		// perform measures on homogeneous configuration and print swaps evolution
		if ( (conf[0].update_index) % param.d_measevery == 0)
		{
			perform_measures_localobs(&(conf[0]), &geo, &param, datafilep, topofilep, &aux_conf);
			print_replicas_labels(swaptrackfilep, conf, &param);
		}

		// save current configurations and backup copies
		if ( param.d_saveconf_backup_every != 0 )
		{
			if ( (conf[0].update_index) % param.d_saveconf_backup_every == 0)
			{
				write_replicas(conf, &param);
				write_replicas_backup(conf, &param);
				write_rng_state(&rng_state, &param);
			}
		}
	}

	// Monte Carlo ends
	time(&finish_date);
	finish_time=clock();
	fprintf(stdout, "#Simulation time: %.10lf s\n", ((double)(finish_time-start_time))/CLOCKS_PER_SEC );

	// save last configurations
	write_replicas(conf, &param);
	write_rng_state(&rng_state, &param);

	// write simulations details on file
	print_simulation_details_cpn(input_file_name, &param, &start_date, &finish_date, start_time, finish_time);

	// close data file
	fclose(datafilep);

	// close topo file
	fclose(topofilep);

	// close swap tracking file
	end_swap_track_file(&swaptrackfilep, &param);

	// print swap acceptances of parallel tempering
	print_swap_acceptances(&swap_counter, &param);

	// free CPN replicas confs
	free_CPN_replicas(conf, &param);

	// free CPN aux conf
	free_CPN_conf(&aux_conf, &param);

	// free swap acceptances
	free_swap_acceptances(&swap_counter, &param);

	// free geometry
	free_geometry(&geo, &param);

	// free rectangles
	free_rectangles_hierarchic_upd(most_update, &param);

	// free rectangle params
	free_param(&param);
}

int main (int argc, char **argv)
{
	char input_file_name[STD_STRING_LENGTH];
	if(argc != 2)
	{
		printf("\n");
		printf("__________________________________________________________________________________________________________________________________\n");
		printf("|________________________________________________________________________________________________________________________________|\n");
		printf("||                                                                                                                              ||\n");
		printf("||                                                                                                                              ||\n");
		printf("||      ,o888888o.    8 888888888o   b.             8      8888888 8888888888 ,o888888o.     8 888888888o       ,o888888o.      ||\n");
		printf("||     8888     `88.  8 8888    `88. 888o.          8            8 8888    . 8888     `88.   8 8888    `88.  . 8888     `88.    ||\n");
		printf("||  ,8 8888       `8. 8 8888     `88 Y88888o.       8            8 8888   ,8 8888       `8b  8 8888     `88 ,8 8888       `8b   ||\n");
		printf("||  88 8888           8 8888     ,88 .`Y888888o.    8            8 8888   88 8888        `8b 8 8888     ,88 88 8888        `8b  ||\n");
		printf("||  88 8888           8 8888.   ,88' 8o. `Y888888o. 8            8 8888   88 8888         88 8 8888.   ,88' 88 8888         88  ||\n");
		printf("||  88 8888           8 888888888P'  8`Y8o. `Y88888o8            8 8888   88 8888         88 8 888888888P'  88 8888         88  ||\n");
		printf("||  88 8888           8 8888         8   `Y8o. `Y8888            8 8888   88 8888        ,8P 8 8888         88 8888        ,8P  ||\n");
		printf("||  `8 8888       .8' 8 8888         8      `Y8o. `Y8            8 8888   `8 8888       ,8P  8 8888         `8 8888       ,8P   ||\n");
		printf("||     8888     ,88'  8 8888         8         `Y8o.`            8 8888    ` 8888     ,88'   8 8888          ` 8888     ,88'    ||\n");
		printf("||      `8888888P'    8 8888         8            `Yo            8 8888       `8888888P'     8 8888             `8888888P'      ||\n");
		printf("||                                                                                                                              ||\n");
		printf("||______________________________________________________________________________________________________________________________||\n");
		printf("|________________________________________________________________________________________________________________________________|\n");
		printf("\n");
		printf("Package: %s-v%s\n", PACKAGE_NAME, PACKAGE_VERSION);
		printf("Compiled from main %s\n", __FILE__);
		printf("Description: lattice simulations of 2d CP^{N-1} models topology via parallel tempering\n\n");
		printf("Author: Claudio Bonanno\n");
		printf("Other contributors: Mario Berni, Davide Vadacchino\n");
		printf("Bug report: %s\n", PACKAGE_BUGREPORT);
		printf("\nCompiled with N = %i\n", N);
		#ifdef __INTEL_COMPILER
		printf("Compiled with icc\n");
		#elif defined( __GNUC__ )
		printf("Compiled with gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		#endif
		printf("Usage: %s input_file \n", argv[0]);
		return(EXIT_FAILURE);
	}
	else
	{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
		{
			fprintf(stderr, "Input file name too long. Increase STD_STRING_LENGTH in include/macro.h\n");
			return(EXIT_FAILURE);
		}
		else
		{
			strcpy(input_file_name, argv[1]);
			real_main(input_file_name);
			return(EXIT_SUCCESS);
		}
	}
}

#endif
