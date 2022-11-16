#ifndef CPN_CONF_CHECK_C
#define CPN_CONF_CHECK_C

#include "../include/macro.h" // NOTE: must be included before stdlib.h to activate posix_memalign correctly
#include "../include/cpn_conf.h"
#include "../include/cpn_param.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void read_CPN_conf_header(CPN_Param *, char const * const);

void real_main(char const * const conf_file_name) 
{
	CPN_Conf conf;
	CPN_Param param;
	
	// get lattice volume from conf header
	read_CPN_conf_header(&param, conf_file_name); // function defined in this file
	// allocate conf
	allocate_CPN_conf(&conf, &param); // function defined in lib/cpn_conf_def.c
	// read conf, compute hash of read conf and compare it with stored hash
	read_CPN_conf_from_file(&conf, &param, conf_file_name); // function defined in lib/cpn_conf_def.c
}

void read_CPN_conf_header(CPN_Param *param, char const * const conf_file_name)
{
	long stored_conf_update_index;
	int stored_size[2];
	int stored_N, err, mu;
	char hash_old[2*MD5_DIGEST_LENGTH+1];
	FILE *fp;

	// open in normal mode and read header
	fp=fopen(conf_file_name, "r"); // open conf file with name conf_file_name
	if(fp==NULL)
	{
		fprintf(stderr, "Error in opening the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	else
	{
		err = fscanf(fp, "%ld %d %d %d", &stored_conf_update_index, &stored_N, &(stored_size[0]), &(stored_size[1]));
		if (err != 4)
		{
			fprintf(stderr, "Error in reading the file %s (%s, %d)\n", conf_file_name, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		else if (stored_N != N)
		{
			fprintf(stderr, "Error: stored value of N %d does not match expected one %d (%s, %d)\n", stored_N, N, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		else
		{
			// store lattice sizes and volume
			param->d_volume=1;
			for (mu=0; mu<2; mu++)
			{
				param->d_size[mu] = stored_size[mu];
				param->d_volume *= stored_size[mu];
			}
		}
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
}

int main (int argc, char **argv)
{
	char conf_file_name[STD_STRING_LENGTH];
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
		printf("Description: tool to check integrity of a 2d CP^{N-1} models field configuration\n\n");
		printf("Author: Claudio Bonanno\n");
		printf("Bug report: %s\n", PACKAGE_BUGREPORT);
		printf("\nCompiled with N = %i\n", N);
		#ifdef __INTEL_COMPILER
			printf("Compiled with icc\n");
		#elif defined( __GNUC__ )
			printf("Compiled with gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
		#endif
		printf("Usage: %s conf_file \n", argv[0]);
		return(EXIT_FAILURE);
	}
	else
	{
		if(strlen(argv[1]) >= STD_STRING_LENGTH)
		{
			fprintf(stderr, "Conf file name too long. Increase STD_STRING_LENGTH in include/macro.h\n");
			return(EXIT_FAILURE);
		}
		else
		{
			#ifdef ENABLE_MD5_HASH
			strcpy(conf_file_name, argv[1]);
			real_main(conf_file_name);
			return(EXIT_SUCCESS);
			#else
			fprintf(stderr, "Error: openssl/md5.h header not found, cannot compute hash of conf %s\n", argv[1]);
			return(EXIT_FAILURE);
			#endif
		}
	}
}

#endif
