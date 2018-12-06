/**
* @file aks_main.c
* Main program file to calculate the AKS algorithm for a given number
* @author Alexandru Paler and Fabio Campos
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>
#include <libspe2.h>
#include <pthread.h>
#include <libmisc.h>
#include <signal.h>
#include <math.h>
#include <malloc_align.h>
#include "gmp.h"
#include "newton_gmp.h"
#include "findR.h"
#include "euler.h"

#define MAX_SPU_THREADS 3

struct timeval timev1, timev2; /* timing-var. */

/*the name of the spu program*/
extern spe_program_handle_t spu_polydiv;

typedef struct ppu_thread_data {
	spe_context_ptr_t ctx;
	void* argp;
} ppu_thread_data_t;

/**
 * Function that is loading the used to launch a spu context
 * @param	pointer to the parameter that is passed to the spu
 * @return 	void
*/
void *ppu_pthread_function(void* arg) {
	unsigned int entry = SPE_DEFAULT_ENTRY;
	ppu_thread_data_t* ptr = (ppu_thread_data_t*)arg;

	if (spe_context_run(ptr->ctx, &entry, 0, ptr->argp, NULL, NULL) < 0) {
		perror("Failed running context");
		exit(1);
	}

	pthread_exit(NULL);
}

/**
 * Function used to launch the work on all the available or needed spu's
 * @param	the number to be tested
 * @param	the R number
 * @param	the limit to test until
 * @param	the number of spus to launch
 * @return 	void
*/
void startSPEs(unsigned int* myInput, unsigned int myRint, unsigned int myLimit, int spuNum) {
	unsigned int message = 0;
	unsigned int msgNull = 0;
	unsigned int isPrime = 1;

	int modCounter;
	int exitCounter=0;
	int result;
	int i;
	int spu_threads;
	int mbox_status;

	ppu_thread_data_t data[spuNum];
	control_block cb[spuNum] __attribute__ ((aligned (128)));

	pthread_t threads[spuNum];

	/*Determine the number of SPE threads to create*/
	spu_threads = spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1);

	if (spu_threads > spuNum) {
		spu_threads = spuNum;
	}
	modCounter = myLimit % spuNum;

	/*Create threads that contain the step3 context*/
	for (i = 0; i < spu_threads; i++) {
		/*Create context*/
		if ((data[i].ctx = spe_context_create(SPE_MAP_PS, NULL)) == NULL) {
			perror("Failed creating context");
			exit(1);
		}
		if (spe_program_load(data[i].ctx, &spu_polydiv)) {
			perror("Failed loading program");
			exit(1);
		}

		/*load the needed information in the structure that is being passed to the context*/
		data[i].argp = &cb[i];
		cb[i].rank = i;
		cb[i].mpm_num[0] = myInput[0];
		cb[i].mpm_num[1] = myInput[1];
		cb[i].mpm_num[2] = myInput[2];
		cb[i].mpm_num[3] = myInput[3];

		cb[i].r = myRint;
		if(modCounter==0) {
			message = message + 1;
			cb[i].from = message;

			message = message + ((myLimit/spuNum)) - 1;
			cb[i].until = message;
		} else {
			message = message + 1;
			cb[i].from = message;

			message = message + (myLimit/spuNum) - 1 + 1;
			cb[i].until = message;
			modCounter = modCounter - 1;
		}

		if ((cb[i].sig1 = spe_ps_area_get(data[i].ctx, SPE_SIG_NOTIFY_1_AREA))
				== NULL) {
			printf("Failed call to spe_ps_area_get %d", i);
			exit(1);
		}

		if ((cb[i].ls = spe_ls_area_get(data[i].ctx)) == NULL) {
			perror("Failed call to spe_getls %d");
			exit(1);
		}

		/*Create and run thread for each SPE context*/
		if (pthread_create(&threads[i], NULL, &ppu_pthread_function, &data[i])) {
			perror("Failed creating thread");
			exit(1);
		}
	}

	/*have all the threads sent their exit message?*/
	while(exitCounter < spuNum) {
		for (i=0; i<spu_threads; i++) {
			if(spe_out_mbox_status(data[i].ctx)>0) {
				spe_out_mbox_read(data[i].ctx, &message, 1);
				if(message == EXIT_MESSAGE) {
					exitCounter++;
				}
				if(message == 0) {
					isPrime = 0;
				}
			}
			message = 1;
		}
	}

	/*Wait for SPU-thread to complete execution*/
	for (i = 0; i < spu_threads; i++) {
		if (pthread_join(threads[i], NULL)) {
			perror("PPU: Failed pthread_join");
			exit(1);
		}
	}

	/*what is the final result?*/
	if(isPrime == 0) {
		printf(", false");
	} else {
		printf(", true");
	}
}

/**
 * Main entry point of the program
 * @param	number of argyuments
 * @return 	0 if everyting ok
*/
int main(int argc, char* argv[]) {

	mpz_t a;
	mpz_t myR;
	int myRint;
	int myEulerCounter;
	int isPrime;
	int logn;
	int myLimit;
	int spuNum = 1;

	isPrime = 1;

	mpz_init(a);
	mpz_init(myR);

	if(argc > 2) {
		if (mpz_set_str(a, argv[1], 0) != 0) {
			/* this is a prototype without error handling */
		}
		spuNum = strtol(argv[2], NULL, 10);
	}

	gmp_printf("%Zd",a);

	/*
		STEP 1 - Is the number a perfect power?
	*/
	gettimeofday(&timev1, NULL);

	if ((newton_it(a))) {
		isPrime = 0;
	}
	gettimeofday(&timev2, NULL);

	/*step1 time*/
	(void) printf(",%d,%f", spuNum, (timev2.tv_sec - timev1.tv_sec) + 0.000001 * (timev2.tv_usec - timev1.tv_usec));

	if(isPrime == 0) {
		return 0;
	}
	gettimeofday(&timev1, NULL);

	/*
		STEP 2 - Search for R
	*/
	myRint = findR(&myR, a);
	(void) printf(",%d", myRint);
	fflush(stdout);
	if(myRint == -1){
		printf(",prime\n");
		return 0;
	}

	gettimeofday(&timev2, NULL);

	/*step2 time*/
	(void) printf(",%f", (timev2.tv_sec - timev1.tv_sec) + 0.000001 * (timev2.tv_usec - timev1.tv_usec));

	/*
		STEP 3 - Intensive Polynom Computations
	*/
	gettimeofday(&timev1, NULL);

	myEulerCounter = euler_probDiv(myRint);
	logn = (int)mpz_sizeinbase(a, 10);
	myLimit = sqrt(myEulerCounter) * logn;

	gettimeofday(&timev1, NULL);

	/* Convert GMP to MPM*/
	/*number of elements of the array*/
	int count = (mpz_sizeinbase(a, 2) + 127)/128;

	/*number of bytes to be allocated*/
	unsigned int* mpm_number =  _malloc_align(count/8, 16);
	unsigned int size;

	/*use the gmp export function*/
	mpz_export(mpm_number, &size, 1, 4, 0, 0, a);

	unsigned int ii;
	/*this is a working cheat
	the elements of the array are not returned on the needed position
	that is why they need to be shifted*/
	if(size != 4)
	{
		for(ii=3; ii>=4-size; ii--)
		{
			mpm_number[ii] = mpm_number[(ii+size)%4];
		}
		for(ii=0; ii<4-size; ii++)
		{
			mpm_number[ii] = 0;
		}
	}
	/* End Convert */

	/*start the work*/
	fflush(stdout);
	startSPEs(mpm_number, myRint, myLimit, spuNum);
	gettimeofday(&timev2, NULL);

	/*step3 time*/
	(void) printf(",%f \n", (timev2.tv_sec - timev1.tv_sec) + 0.000001 * (timev2.tv_usec - timev1.tv_usec));

	return 0; /* everything is allright */
}
