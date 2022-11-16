#ifndef CPN_PARAM_C
#define CPN_PARAM_C

#include "../include/cpn_param.h"

void remove_white_lines_and_comments(FILE *input_fp)
{
	int temp_i;
	temp_i=getc(input_fp);
	if(temp_i=='\n' || temp_i==' ' || temp_i=='\043') // scan for white lines and comments
	{
		ungetc(temp_i, input_fp);

		temp_i=getc(input_fp);
		if(temp_i=='\n' || temp_i==' ') // white line
		{
			do
			{
				temp_i=getc(input_fp);
			}
			while(temp_i=='\n' || temp_i==' ');
		}
		ungetc(temp_i, input_fp);

		temp_i=getc(input_fp);
		if(temp_i=='\043')  // comment
		{
			do
			{
				temp_i=getc(input_fp);
			}
			while(temp_i!='\n');
		}
		else
		{
			ungetc(temp_i, input_fp);
		}

		remove_white_lines_and_comments(input_fp);
	}
	else
	{
		ungetc(temp_i, input_fp);
	}
}

void read_input(char const * const input_file_name, CPN_Param *param)
{
	FILE *input_fp;
	char str[STD_STRING_LENGTH], temp_str[STD_STRING_LENGTH];
	double temp_d;
	int temp_i, i, mu, err, end=1;
	long temp_l;

	// to avoid problems when computing param->d_ngrid
	param->d_grid_step=1;
	param->d_grid_max=0;

	input_fp=fopen(input_file_name, "r");  // open the input file
	if(input_fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		while(end==1) // slide the file
		{
			remove_white_lines_and_comments(input_fp);
			err=fscanf(input_fp, "%s", str);
			if(err!=1)
			{
				fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}

			if(strncmp(str, "size", 4)==0)
			{
				for(i=0; i<2; i++)
				{
					err=fscanf(input_fp, "%d", &temp_i);
					if(err!=1)
					{
						fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
						exit(EXIT_FAILURE);
					}
					param->d_size[i]=temp_i;
				}
			}

			else if(strncmp(str, "beta", 4)==0)
			{ 
				err=fscanf(input_fp, "%lf", &temp_d);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_beta=temp_d;
			}
			else if(strncmp(str, "theta", 5)==0)
			{
				err=fscanf(input_fp, "%lf", &temp_d);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_theta=temp_d;
			}
			else if(strncmp(str, "grid_step", 9)==0)
			{
				err=fscanf(input_fp, "%lf", &temp_d);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_grid_step=temp_d;
			}
			else if(strncmp(str, "grid_max", 8)==0)
			{
				err=fscanf(input_fp, "%lf", &temp_d);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_grid_max=temp_d;
			}
			else if(strncmp(str, "num_single_site_stoc_upd", 24)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_stoc_single_site_upd=temp_i;
			}
			else if(strncmp(str, "num_MC_step", 11)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_MC_step=temp_i;
			}
			else if(strncmp(str, "meas_every", 10)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_measevery=temp_i;
			}
			else if(strncmp(str, "num_micro", 9)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_num_micro=temp_i;
			}
			else if(strncmp(str, "num_norm", 8)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_num_norm=temp_i;
			}
			else if(strncmp(str, "start", 5)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_start=temp_i;
			}
			else if(strncmp(str, "save_conf_every", 15)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_saveconf_backup_every=temp_i;
			}
			else if(strncmp(str, "num_cool_step", 13)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_coolsteps=temp_i;
			}
			else if(strncmp(str, "meas_cool_every", 15)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_coolevery=temp_i;
			}
			else if(strncmp(str, "conf_file", 9)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_conf_file, temp_str);
			}
			else if(strncmp(str, "data_file", 9)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_data_file, temp_str);
			}
			else if(strncmp(str, "topo_file", 9)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_topo_file, temp_str);
			}
			else if(strncmp(str, "log_file", 8)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_log_file, temp_str);
			}
			else if(strncmp(str, "rng_seed", 8)==0)
			{ 
				err=fscanf(input_fp, "%ld", &temp_l);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_seed=temp_l;
			}	
			else if(strncmp(str, "rng_start", 9)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_rng_start=temp_i;
			}
			else if(strncmp(str, "rng_state_file", 14)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_rng_file, temp_str);
			}
			else if(strncmp(str, "defect_size", 11)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_L_defect=temp_i;
			}
			else if(strncmp(str, "num_replica", 11)==0)
			{
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_N_replica_pt=temp_i;
			}
			else if(strncmp(str, "swap_acc_file", 13)==0)
			{ 
			err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_swap_accept_file, temp_str);
			}  
			else if(strncmp(str, "swap_track_file", 15)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_swap_tracking_file, temp_str);
			}
			else if(strncmp(str, "topo_potential_file", 19)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_topo_potential_file, temp_str);
			}     
			else if(strncmp(str, "multicanonic_acc_file", 21)==0)
			{ 
				err=fscanf(input_fp, "%s", temp_str);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				strcpy(param->d_multicanonic_acc_file, temp_str);
			}   
			else if(strncmp(str, "hierarc_upd", 11)==0)
			{ 
				err=fscanf(input_fp, "%d", &temp_i);
				if(err!=1)
				{
					fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
					exit(EXIT_FAILURE);
				}
				param->d_N_hierarc_levels=temp_i;
				  
				if (temp_i > 0) // if performing hierarchic updates, read the other parameters
				{
					err=posix_memalign( (void **) &(param->d_L_rect), (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(int));
					if (err!=0)
					{
						fprintf(stderr, "Problems in allocating hierarchical update parameters! (%s, %d)\n", __FILE__, __LINE__);
						exit(EXIT_FAILURE);
					}
					for(i=0;i<param->d_N_hierarc_levels;i++)
					{
						err=fscanf(input_fp, "%d", &temp_i);
						if(err!=1)
						{
							fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
							exit(EXIT_FAILURE);
						}				
						param->d_L_rect[i]=temp_i;
					}
									
					err=posix_memalign( (void **) &(param->d_N_sweep_rect), (size_t) INT_ALIGN, (size_t) param->d_N_hierarc_levels * sizeof(int));
					if (err!=0)
					{
						fprintf(stderr, "Problems in allocating hierarchical update parameters! (%s, %d)\n", __FILE__, __LINE__);
						exit(EXIT_FAILURE);
					}
					for(i=0;i<param->d_N_hierarc_levels;i++)
					{
						err=fscanf(input_fp, "%d", &temp_i);
						if(err!=1)
						{
							fprintf(stderr, "Error in reading the file %s (%s, %d)\n", input_file_name, __FILE__, __LINE__);
							exit(EXIT_FAILURE);
						}
						param->d_N_sweep_rect[i]=temp_i;
					}
				}
			}
			else
			{
				fprintf(stderr, "Error: unrecognized option %s in the file %s (%s, %d)\n", str, input_file_name, __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}

			remove_white_lines_and_comments(input_fp);

			// check if the read line is the last one
			temp_i=getc(input_fp);
			if (temp_i==EOF) end=0; // if end of file is reached, end=0 and exit the while loop
			else ungetc(temp_i, input_fp);
		} // closes initial while

		fclose(input_fp);

		// VARIOUS CHECKS
		if(param->d_L_defect>param->d_size[1]) // defect is along the space direction
		{
			fprintf(stderr, "Error: defect_size is > size[1], must be smaller (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		
		if(param->d_N_replica_pt<1)
		{
			fprintf(stderr, "Error: num_replica parameter must be > 0 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
				
		if( (param->d_start!=0) && (param->d_start!=1) && (param->d_start!=2) )
		{
			fprintf(stderr, "Error: start parameter must be either 0, 1 or 2 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		if (param->d_measevery<1)
		{
			fprintf(stderr, "Error: meas_every parameter must be > 0 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);		
		}

		if (param->d_num_norm<1)
		{
			fprintf(stderr, "Error: num_norm parameter must be > 0 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);		
		}

		if (param->d_seed==0)
		{
			fprintf(stderr, "Error: rng_seed parameter must be != 0 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		if ( (param->d_rng_start!=0) && (param->d_rng_start!=1) )
		{
			fprintf(stderr, "Error: rng_start parameter must be either 0 or 1 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	
		err=0;
		for(mu=0; mu<2; mu++)
		{
			if(param->d_size[mu]==1) err=1;
		}
		if(err==1)
		{
			fprintf(stderr, "Error: all sizes have to be larger than 1 (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}

		init_derived_constants(param);
	} // closes initial if
}

// compute derived constants	
void init_derived_constants(CPN_Param *param)
{
	int mu;
	param->d_n_grid=(int)((2.0*param->d_grid_max/param->d_grid_step)+1.0);
	param->d_volume=1;
	for(mu=0; mu<2; mu++) (param->d_volume)*=(param->d_size[mu]);
}

void free_param(CPN_Param *param)
{
	if (param->d_N_hierarc_levels>0)
	{
		free(param->d_L_rect);
		free(param->d_N_sweep_rect);
	}
}
  
// initialize data file
void init_data_file(FILE **dataf, CPN_Param const * const param)
{
	*dataf=fopen(param->d_data_file, "r");
	if(*dataf!=NULL) // file exists
	{
		fclose(*dataf);
		*dataf=fopen(param->d_data_file, "a");
	}
	else // file doesn't exist
	{
		int mu;
		*dataf=fopen(param->d_data_file, "w");
		fprintf(*dataf, "# %f ", param->d_beta);
		for(mu=0; mu<2; mu++) fprintf(*dataf, "%d ", param->d_size[mu]);
		fprintf(*dataf, "\n");
	}
	fflush(*dataf);
}
  
// initialize topo data file
void init_topo_file(FILE **topofilep, CPN_Param const * const param)
{ 
	*topofilep=fopen(param->d_topo_file, "r");
	if(*topofilep!=NULL) // file exists
	{
		fclose(*topofilep);
		*topofilep=fopen(param->d_topo_file, "a");
	}
	else
	{
		int mu;
		*topofilep=fopen(param->d_topo_file, "w");
		fprintf(*topofilep, "# %f ", param->d_beta);
		for(mu=0; mu<2; mu++) fprintf(*topofilep, "%d ", param->d_size[mu]);
		fprintf(*topofilep, "\n");
	}
	fflush(*topofilep);
}

// print simulations details of cpn
void print_simulation_details_cpn(char const * const input_file_name, CPN_Param const * const param, time_t const * const start_date, time_t const * const finish_date,
                                  clock_t const start_time, clock_t const finish_time)
{
	FILE *fp;
	int i;
	
	fp=fopen(param->d_log_file, "w");
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_log_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(fp, "+----------------------------------------------+\n");
		fprintf(fp, "|          Simulation details for cpn          |\n");
		fprintf(fp, "+----------------------------------------------+\n\n");

		fprintf(fp, "Theory: 2d CP^{N-1} model with N = %d\n\n", N);
		fprintf( fp, "Input parameters read from file %s\n\n", input_file_name);
		fprintf(fp, "Lattice size: %d x %d\n", param->d_size[0], param->d_size[1]);
		if (param->d_N_replica_pt>1) fprintf(fp, "Number of replicas used for parallel tempering: %d\n", param->d_N_replica_pt);
		if (param->d_start==0) fprintf(fp, "Simulation started from random cold conf\n");
		if (param->d_start==1) fprintf(fp, "Simulation started from random hot conf\n");
		if (param->d_start==2) fprintf(fp, "Simulation started from stored conf read from file %s", param->d_conf_file);
		if (param->d_N_replica_pt>1) fprintf(fp, "_replica_*\n");
		else fprintf(fp, "\n");

		fprintf(fp, "\n");

		fprintf(fp, "beta:  %.10lf\n", param->d_beta);
		fprintf(fp, "theta: %.10lf\n", param->d_theta);

		if (param->d_L_defect>0) fprintf(fp, "Defect size: %d\n", param->d_L_defect);
		if(param->d_N_hierarc_levels>0)
		{
			fprintf(fp,"Number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
			fprintf(fp,"Rectangles extension: ");
			for(i=0;i<param->d_N_hierarc_levels;i++) fprintf(fp,"%d ", param->d_L_rect[i]);
			fprintf(fp,"\n");
			fprintf(fp,"Number of sweeps per hierarchic level: ");
			for(i=0;i<param->d_N_hierarc_levels;i++) fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
		}
		fprintf(fp, "\n");
	
		if (param->d_N_replica_pt==1) fprintf(fp, "Number of standard local updating steps: %d\n", param->d_MC_step);
		else fprintf(fp, "Number of parallel tempering steps: %d\n", param->d_MC_step);
		fprintf(fp, "Over-relaxation/over-heat-bath ratio: %d / 1\n", param->d_num_micro);
		fprintf(fp, "Measures taken every: %d\n", param->d_measevery);
		fprintf(fp, "Conf normalized every: %d\n", param->d_num_norm);
		fprintf(fp, "Conf saved every: %d\n", param->d_saveconf_backup_every);
		fprintf(fp, "\n");

		fprintf(fp, "Max num of cooling steps: %d\n", param->d_coolsteps);
		fprintf(fp, "Topo obs measured every %d cooling steps\n", param->d_coolevery);
		fprintf(fp, "\n");

		if (param->d_rng_start==0) fprintf(fp, "Rng initializated from seed %ld\n", param->d_seed);
		if (param->d_rng_start==1) fprintf(fp, "Rng initializated from last rng state read from file %s (seed = %ld)\n", param->d_rng_file, param->d_seed);
		fprintf(fp, "\n");
		
		if (is_little_endian() == 0) fprintf(fp, "Little endian machine\n");
		else fprintf(fp, "Big endian machine\n");
		fprintf(fp, "\n");

		fprintf(fp, "Simulation start: %s", ctime(start_date));
		fprintf(fp, "Simulation end:   %s", ctime(finish_date));
		fprintf(fp, "Simulation time:  %.5lf seconds\n", ((double)(finish_time-start_time))/CLOCKS_PER_SEC );
		fprintf(fp, "\n");
	
		fclose(fp);
	}
}

// print simulations details of multicanonic_cpn
void print_simulation_details_multicanonic_cpn(char const * const input_file_name, CPN_Param const * const param, time_t const * const start_date, time_t const * const finish_date,
                                               clock_t const start_time, clock_t const finish_time)
{
	FILE *fp;
	int i;
	
	fp=fopen(param->d_log_file, "w");
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_log_file, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(fp, "+----------------------------------------------+\n");
		fprintf(fp, "|   Simulation details for multicanoninc_cpn   |\n");
		fprintf(fp, "+----------------------------------------------+\n\n");

		fprintf(fp, "Theory: 2d CP^{N-1} model with N = %d\n\n", N);
		fprintf( fp, "Input parameters read from file %s\n\n", input_file_name);
		fprintf(fp, "Lattice size: %d x %d\n", param->d_size[0], param->d_size[1]);
		if (param->d_N_replica_pt>1) fprintf(fp, "Number of replicas used for parallel tempering: %d\n", param->d_N_replica_pt);
		if (param->d_start==0) fprintf(fp, "Simulation started from random cold conf\n");
		if (param->d_start==1) fprintf(fp, "Simulation started from random hot conf\n");
		if (param->d_start==2) fprintf(fp, "Simulation started from stored conf read from file %s", param->d_conf_file);
		if (param->d_N_replica_pt>1) fprintf(fp, "_replica_*\n");
		else fprintf(fp, "\n");

		fprintf(fp, "\n");

		fprintf(fp, "beta:  %.10lf\n", param->d_beta);
		fprintf(fp, "theta: %.10lf\n", param->d_theta);

		if (param->d_L_defect>0) fprintf(fp, "Defect size: %d\n", param->d_L_defect);
		if(param->d_N_hierarc_levels>0)
		{
			fprintf(fp,"Number of hierarchical levels: %d\n", param->d_N_hierarc_levels);
			fprintf(fp,"Rectangles extension: ");
			for(i=0;i<param->d_N_hierarc_levels;i++) fprintf(fp,"%d ", param->d_L_rect[i]);
			fprintf(fp,"\n");
			fprintf(fp,"Number of sweeps per hierarchic level: ");
			for(i=0;i<param->d_N_hierarc_levels;i++) fprintf(fp,"%d ", param->d_N_sweep_rect[i]);
		}
		fprintf(fp, "\n\n");
	
		if (param->d_N_replica_pt==1) fprintf(fp, "Number of standard local stochastic updating steps: %d\n", param->d_MC_step);
		else fprintf(fp, "Number of multicanonic stochastic parallel tempering steps:    %d\n", param->d_MC_step);
		fprintf(fp, "Number of multicanonic stochastic single site/link refreshing steps for each updating step: %d\n", param->d_stoc_single_site_upd);
		fprintf(fp, "Topo-potential read from file %s\n", param->d_topo_potential_file);
		fprintf(fp, "Topo-potential defined on the interaval [-%lf, %lf] in steps of %lf\n", param->d_grid_max, param->d_grid_max, param->d_grid_step);
		fprintf(fp, "Over-relaxation/over-heat-bath ratio: %d / 1\n", param->d_num_micro);
		fprintf(fp, "Measures taken every: %d\n", param->d_measevery);
		fprintf(fp, "Conf normalized every: %d\n", param->d_num_norm);
		fprintf(fp, "Conf saved every: %d\n", param->d_saveconf_backup_every);
		fprintf(fp, "\n");

		fprintf(fp, "Max num of cooling steps: %d\n", param->d_coolsteps);
		fprintf(fp, "Topo obs measured every %d cooling steps\n", param->d_coolevery);
		fprintf(fp, "\n");

		if (param->d_rng_start==0) fprintf(fp, "Rng initializated from seed %ld\n", param->d_seed);
		if (param->d_rng_start==1) fprintf(fp, "Rng initializated from last rng state read from file %s (seed = %ld)\n", param->d_rng_file, param->d_seed);
		fprintf(fp, "\n");
		
		if (is_little_endian() == 0) fprintf(fp, "Little endian machine\n");
		else fprintf(fp, "Big endian machine\n");
		fprintf(fp, "\n");

		fprintf(fp, "Simulation start: %s", ctime(start_date));
		fprintf(fp, "Simulation end:   %s", ctime(finish_date));
		fprintf(fp, "Simulation time:  %.5lf seconds\n", ((double)(finish_time-start_time))/CLOCKS_PER_SEC );
		fprintf(fp, "\n");
	
		fclose(fp);
	}
}

#endif
