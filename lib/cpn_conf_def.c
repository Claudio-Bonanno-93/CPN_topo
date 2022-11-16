#ifndef CPN_CONF_DEF_C
#define CPN_CONF_DEF_C

#include "../include/cpn_conf.h"
  
// allocate all replicas for parallel tempering
void init_CPN_replicas(CPN_Conf **conf, CPN_Param const * const param, RNG_Param *rng_state)
{
	int i=0, err;
	char conf_file_name[STD_STRING_LENGTH], r[STD_STRING_LENGTH];
	strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file

	// allocate the vector to store replicas
	err=posix_memalign((void **) conf, (size_t) DOUBLE_ALIGN, (size_t) param->d_N_replica_pt * sizeof(CPN_Conf));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the parallel tempering replicas! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	if (param->d_N_replica_pt==1)
	{
		allocate_CPN_conf(&((*conf)[i]),param); // allocate memory to store CPN conf
		init_CPN_conf(&((*conf)[i]), param, conf_file_name, rng_state); // initialize CPN conf
		init_bound_cond(&((*conf)[i]),i,param); // initialize boundary conditions parameters
		((*conf)[i]).conf_label=i;
	}
	else
	{
		for(i=0; i<param->d_N_replica_pt; i++)
		{
			sprintf(r, "%d", i); // r='i'
			strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file
			strcat(conf_file_name, "_replica_"); // conf_file_name = param->d_conf_file + '_replica_'
			strcat(conf_file_name, r); // conf_file_name = param->d_conf_file + '_replica_i'

			allocate_CPN_conf(&((*conf)[i]), param);
			init_CPN_conf(&((*conf)[i]), param, conf_file_name, rng_state);
			init_bound_cond(&((*conf)[i]), i, param);
			((*conf)[i]).conf_label=i;
		}
	}
}

// allocate memory for a CPN conf
void allocate_CPN_conf(CPN_Conf *conf, CPN_Param const * const param)
{
	long i;
	int err;
	err=posix_memalign((void**) &(conf->z), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(cmplx *));
	if(err!=0)
    {
		fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
    }
	err=posix_memalign((void**) &(conf->U), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(cmplx *));
	if(err!=0)
    {
		fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
    }
	for(i=0; i<param->d_volume; i++)
	{
		err=posix_memalign((void**) &(conf->z[i]), (size_t) DOUBLE_ALIGN, (size_t) N * sizeof(cmplx));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		err=posix_memalign((void**) &(conf->U[i]), (size_t) DOUBLE_ALIGN, (size_t) 2 * sizeof(cmplx));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}
}

// initialize single CPN conf replica
void init_CPN_conf(CPN_Conf *conf, CPN_Param const * const param, char const * const conf_file_name, RNG_Param *rng_state) 
{ 
	int i, mu;
	// random cold (ordered) initialization
	if (param->d_start == 0)
	{
		conf->update_index = 0;
		for ( i=0 ; i<param->d_volume ; i++ ) 
		{
			for( mu=0 ; mu<2 ; mu++) conf->U[i][mu] = 1.0 + I * 0.0;
			vector_random_cold(conf->z[i], rng_state);
		}
	}

	// random hot (casual) initialization
	if (param->d_start == 1)
	{
		conf->update_index = 0;
		for (i=0; i<param->d_volume; i++) 
		{
			for(mu=0; mu<2; mu++)
			{
				conf->U[i][mu] = cmplx_rand_num(rng_state);
				conf->U[i][mu] /= cmplx_abs( conf->U[i][mu]);
			}
			vector_random_hot(conf->z[i], rng_state);
		}
	}
	
	// read conf from file
	if (param->d_start == 2) read_CPN_conf_from_file(conf, param, conf_file_name);
}

// initialization of the defect for the replica with label a
void init_bound_cond(CPN_Conf *conf, int const a, CPN_Param const * const param)
{
	int err, mu;
	long i;

	//allocation of C[i][mu]
	err=posix_memalign((void**) &(conf->C), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double *));
	if(err!=0)
	{
		fprintf(stderr, "Problems in allocating the defect!\n");
		exit(EXIT_FAILURE);
	}

	for(i=0; i<param->d_volume; i++)
	{
		err=posix_memalign((void**) &(conf->C[i]), (size_t) DOUBLE_ALIGN, (size_t) 2 * sizeof(double));
		if(err!=0)
		{
			fprintf(stderr, "Problems in allocating the defect!\n");
			exit(EXIT_FAILURE);
		}
	}

	//initialization of C[r][mu]
	for(i=0; i<param->d_volume; i++)
	{
		for(mu=0; mu<2; mu++)
		{
			conf->C[i][mu]=1.0;
			// defect affects link along the mu=0 direction being set along the mu=1 direction
			// if mu=0, modify boundary conditions on the defect if more than 1 replica is used
			if ( (param->d_N_replica_pt>1) && (mu==0) && (is_on_defect(i, param) == 0) )
				conf->C[i][mu] -= (((double)(a))/(((double)(param->d_N_replica_pt-1)))); // C = 1 - a / (N_replica - 1)
		}
	}
}

// check if a site lies on the defect: 0 = is on defect, 1 = is not on defect
int is_on_defect(long const i, CPN_Param const * const param)
{
	long cart_coord[2];
	si_to_cart(cart_coord, i, param);
	if( (cart_coord[0]==(param->d_size[0]-1)) && (cart_coord[1]<param->d_L_defect) ) return 0;
	else return 1;
}

// normalize all parallel tempering replicas
void normalize_replicas(CPN_Conf *conf, CPN_Param const * const param)
{
	int i;
	for (i=0; i<param->d_N_replica_pt; i++) normalize_CPN_conf(&(conf[i]),param);
}

// normalize a single configurations
void normalize_CPN_conf(CPN_Conf *conf, CPN_Param const * const param) 
{ 
	int mu; 
	long i;

	for (i=0; i<param->d_volume; i++) 
	{
		vector_normalization(conf->z[i]);
		for(mu=0; mu<2; mu++) conf->U[i][mu] /= cmplx_abs(conf->U[i][mu]);
	}
} 

// copy conf in auxiliary conf
void copyconf(CPN_Conf const * const conf, CPN_Param const * const param, CPN_Conf * aux_conf) 
{ 
	int i, mu;
	for (i=0; i<param->d_volume; i++)
	{ 
		vector_equal(aux_conf->z[i], conf->z[i]); // aux_conf->z = conf->z
		for (mu=0; mu<2; mu++) aux_conf->U[i][mu] = conf->U[i][mu]; // aux_conf->U = conf->U
	} 
} 

// save replicas confs
void write_replicas(CPN_Conf const * const conf, CPN_Param const * const param)
{
	int i=0;
	char conf_file_name[STD_STRING_LENGTH], r[STD_STRING_LENGTH];
	strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file

	if(param->d_N_replica_pt==1) write_CPN_conf_on_file(&(conf[i]), param, conf_file_name);
	else
	{
		for(i=0; i<param->d_N_replica_pt; i++)
		{
			sprintf(r, "%d", i); // r='i'
			strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file
			strcat(conf_file_name, "_replica_"); // conf_file_name = param->d_conf_file + '_replica_'
			strcat(conf_file_name, r); // conf_file_name = param->d_conf_file + '_replica_i'
			write_CPN_conf_on_file(&(conf[i]), param, conf_file_name);
		}
	}
}

// save replicas confs for backup
void write_replicas_backup(CPN_Conf const * const conf, CPN_Param const * const param)
{
	static int counter=0;
	int i=0;
	char conf_file_name[STD_STRING_LENGTH], r[STD_STRING_LENGTH], backup_count[STD_STRING_LENGTH];
	sprintf(backup_count, "_back_%d", counter);
	strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file
	strcat(conf_file_name, backup_count); // conf_file_name = param->d_conf_file + '_back_counter'
	
	if(param->d_N_replica_pt==1) write_CPN_conf_on_file(&(conf[i]), param, conf_file_name);
	else
	{
		for(i=0; i<param->d_N_replica_pt; i++)
		{
			sprintf(r, "%d", i); // r='i'
			strcpy(conf_file_name, param->d_conf_file); // conf_file_name = param->d_conf_file
			strcat(conf_file_name, "_replica_"); // conf_file_name = param->d_conf_file + '_replica_'
			strcat(conf_file_name, r); // conf_file_name = param->d_conf_file + '_replica_i'
			strcat(conf_file_name, backup_count); // conf_file_name = param->d_conf_file + '_replica_i' + '_back_counter'
			write_CPN_conf_on_file(&(conf[i]), param, conf_file_name);
		}
	}
	counter = 1 - counter; // alternate between "_back_0" and "_back_1" in order to keep two different backups of the configurations
}

// save conf on file with name
void write_CPN_conf_on_file(CPN_Conf const * const conf, CPN_Param const * const param, char const * const conf_file_name)
{
	long i;
	int mu;
	FILE *fp;
	#ifdef ENABLE_MD5_HASH
	char hash[2*MD5_DIGEST_LENGTH+1];
	#endif

	// open in normal mode and write header (first line)
	fp=fopen(conf_file_name, "w"); // open conf file with name conf_file_name
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
    {
		fprintf(fp, "%ld %d %d %d", conf->update_index, N, param->d_size[0], param->d_size[1]);
		#ifdef ENABLE_MD5_HASH
		compute_MD5_hash_conf(hash, conf, param);
		fprintf(fp, " %s", hash);
		#endif
		fprintf(fp, "\n");
		fclose(fp);
    }
	
	// open in binary mode and write conf
	fp=fopen(conf_file_name, "ab"); // open conf file with name conf_file_name
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		for (i=0; i<param->d_volume; i++)
		{
			write_cmplx_vec_on_file_binary(fp, conf->z[i]);
			for (mu=0; mu<2; mu++) write_cmplx_num_on_file_binary(fp, &(conf->U[i][mu]));
		}
		fclose(fp);
	}
}

// read conf from file with name
void read_CPN_conf_from_file(CPN_Conf *conf, CPN_Param const * const param, char const * const conf_file_name)
{
	long i;
	int mu, err, stored_N, stored_size_0, stored_size_1;
	long stored_conf_update_index;
	FILE * fp;
	#ifdef ENABLE_MD5_HASH
	char hash_old[2*MD5_DIGEST_LENGTH+1], hash_new[2*MD5_DIGEST_LENGTH+1];
	#endif
	
	// open in normal mode and read header (first line)
	fp=fopen(conf_file_name, "r"); // open conf file with name conf_file_name
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		err = fscanf(fp, "%ld %d %d %d", &stored_conf_update_index, &stored_N, &stored_size_0, &stored_size_1);
		if (err != 4)
		{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		else if( (stored_N != N) || (stored_size_0 != param->d_size[0]) || (stored_size_1 != param->d_size[1]) )
		{
			fprintf(stderr, "Error in reading the file %s: stored parameters (N=%d, L_0=%d, L_1=%d) do not match expected ones (N=%d, L_0=%d, L_1=%d) (%s, %d)\n",
							conf_file_name, stored_N, stored_size_0, stored_size_1, N, param->d_size[0], param->d_size[1], __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		else conf->update_index = stored_conf_update_index;

		#ifdef ENABLE_MD5_HASH
		err = fscanf(fp, "%s", hash_old);
		if (err != 1)
		{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		#endif
	}
	fclose(fp);

	// open in binary mode and read the conf
	fp=fopen(conf_file_name, "rb");
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		// skip header (first line)
		err=0;
		while(err!='\n') err=fgetc(fp);
		// read conf
		for (i=0; i<param->d_volume; i++)
		{
			read_cmplx_vec_from_file_binary(fp, conf->z[i]);
			for(mu=0; mu<2; mu++) read_cmplx_num_from_file_binary(fp, &(conf->U[i][mu]));
		}
		fclose(fp);
	}

	#ifdef ENABLE_MD5_HASH
	compute_MD5_hash_conf(hash_new, conf, param);
	if(strncmp(hash_old, hash_new, (2*MD5_DIGEST_LENGTH+1) )!=0)
	{
		fprintf(stderr, "The computed hash %s does not match the stored one %s (%s, %d)\n", hash_new, hash_old, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	#endif
}

// compute MD5 hash of a given conf
void compute_MD5_hash_conf(char *hash, CPN_Conf const * const conf, CPN_Param const * const param)
{
	#ifdef ENABLE_MD5_HASH
	MD5_CTX mdContext;
	unsigned char c[MD5_DIGEST_LENGTH];
	long i;
	int mu,k;

	// init hash
	MD5_Init(&mdContext);

	// compute hash
	for (i=0; i<param->d_volume; i++)
	{
		for(k=0; k<N; k++) MD5_Update(&mdContext, &(conf->z[i][k]), sizeof(cmplx));
		for(mu=0; mu<2; mu++) MD5_Update(&mdContext, &(conf->U[i][mu]), sizeof(cmplx));
	}

	// finalize hash
	MD5_Final(c, &mdContext);

	// save hash
	for(k=0; k<MD5_DIGEST_LENGTH; k++) sprintf(&(hash[2*k]), "%02x", c[k]);
	#else
	// to avoid compiler warnings about unused variables
	(void) hash;
	(void) conf;
	(void) param;
	#endif
}


void free_CPN_replicas(CPN_Conf *conf, CPN_Param const * const param)
{
	int i;
	for(i=0; i<param->d_N_replica_pt; i++)
	{
		free_CPN_conf(&(conf[i]), param);
		free_bound_cond(&(conf[i]), param);
	}
	free(conf);
}
	
void free_bound_cond(CPN_Conf *conf, CPN_Param const * const param)
{
	long i;
	for(i=0; i<(param->d_volume); i++) free(conf->C[i]);
	free(conf->C);
}
  
void free_CPN_conf(CPN_Conf *conf, CPN_Param const * const param)
{
	long i;
	for(i=0; i<(param->d_volume); i++)
	{
		free(conf->z[i]);
		free(conf->U[i]);
	}
	free(conf->z);
	free(conf->U);
}

#endif
