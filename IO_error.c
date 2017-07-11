#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO_error.h"

extern struct parameter simParams;

void printWarning(char *msg, int level)
{
#if OUTPUT==1
	char str[1024];
	#if PARALLEL==1
		#pragma omp critical
	#endif
	{
		sprintf(str, "out_error/%s_error.log", simParams.fname);
		if( (simParams.fp_errorOut = fopen(str, "a")) == NULL){
			printf("Cannot open error.log\n");
			exit(1);
		}
		/* level allows for different output situations--currently unused */
		fprintf(simParams.fp_errorOut, "time: %d Warning: %s\n",simParams.currentTimestep,msg);
		fclose(simParams.fp_errorOut);
		simParams.numErrors++;
	}
#endif
}
void printError(char *msg, int level)
{
#if OUTPUT==1
	char str[1024];
	#if PARALLEL==1
		#pragma omp critical
	#endif
	{
		sprintf(str, "out_error/%s_error.log", simParams.fname);
		if( (simParams.fp_errorOut = fopen(str, "a")) == NULL){
			printf("Cannot open error.log\n");
			exit(1);
		}
		printf("time: %d Error: %s\n",simParams.currentTimestep,msg);
		fprintf(simParams.fp_errorOut, "time: %d Error: %s\n",simParams.currentTimestep,msg);
		fclose(simParams.fp_errorOut);
		simParams.numErrors++;
	}
#endif
}

void printFatalError(char *msg, int level)
{
	char str[1024];

	#if PARALLEL==1
		#pragma omp critical
	#endif
	{
		sprintf(str, "out_error/%s_error.log", simParams.fname);
		if( (simParams.fp_errorOut = fopen(str, "a")) == NULL){
			printf("Cannot open error.log\n");
			exit(1);
		}
		fprintf(simParams.fp_errorOut, "time: %d Fatal Error: %s\n",simParams.currentTimestep,msg);
		printf("time: %d Fatal Error: %s\n",simParams.currentTimestep,msg);
		fclose(simParams.fp_errorOut);
	}
	//clean up output files
	closeOutputFiles();
	exit(1);
}
