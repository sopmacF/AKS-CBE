/** 
* @file spu_polydiv.c
* The SPU program for the calculations of Step3.
* This program is being loaded on every SPU to compute the following equation 
* (x + a)^n mod (x^r - 1) != x^(n mod r) + a
* @line This has been ported from the corresponding GMP implementation
* @warning The program is able to operate on numbers of maximum 128bits although the implementation has been started to 
* support larger numbers
* The parameters to the class are being sent over using a struct @link signal.h @endlink
* @author Alexandru Paler and Fabio Campos
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <libmpm.h>
#include <mpm_macros.h>
#include <malloc_align.h>
#include <free_align.h>
#include "../signal.h"

#define TRUE_MPM 	0xFFFFFFFF
#define FALSE_MPM 	0x00000000
#define MAX_DEGREE	5000
#define MAX_SIZE	1

typedef struct  {
	int degree __attribute__ ((aligned(16))); /* this value is the real-poly-degree + 1
	 * for example: the poly x^5 has a degree value of 6
	 * in this struct
	 */
	vector unsigned int** coeff; /* Array presenting the coefficients of the polynom */
	int* coeff_size; /* Array containing the sizes of the polynom coefficients */

} poly;

inline void polyPrint(poly* in);

inline void polyExp_Mod_int(poly* resultPoly, int r);

inline void polyMult_Mod(poly* resultPoly, poly* left, poly* right, int r);

inline void polyMult(poly* resultPoly, poly* left, poly* right);

inline int polyIsEqual(poly* left, vector unsigned int* right_coeff1, vector unsigned int* right_coeff0);

inline void polyCopy(poly* ResultPoly, poly* inPoly);

inline int testBitVector(vector unsigned int* number, int number_size, int bit);

control_block cb __attribute__ ((aligned (128)));

vector unsigned int incNr = (vector unsigned int){0, 0, 0, 1};
int incNr_size = 1;

vector unsigned int zeroNr = (vector unsigned int){0, 0, 0, 0};
int zeroNr_size = 1;

vector unsigned int expv;
int expv_size;

vector unsigned int*  temporal;
int temporal_size;

vector unsigned int*  temporal2;
int temporal2_size;

int myLimit;
unsigned int modR;

poly* dummyPoly;
poly* basisPoly;

int rightPoly_degree;
vector unsigned int rightPoly_coeffmodr;
vector unsigned int rightPoly_coeff0;


inline void initVector(poly* poly)
{
	int i __attribute__ ((aligned(16)));
	int j __attribute__ ((aligned(16)));
	for(i=0; __builtin_expect(i < poly->degree, 1); i++)
	{
		poly->coeff_size[i] = zeroNr_size;
		/*for(j=0; j<MAX_SIZE; j++)
		{
			poly->coeff[i][j] = zeroNr;
		}*/
		//poly->coeff[i][0] = zeroNr;
		poly->coeff[i][0] = (vector unsigned int)spu_splats(0);
	}
	poly->degree = 0;
}

inline void initPoly(poly* polyn)
{
	int i __attribute__ ((aligned(16)));

	polyn->coeff = (vector unsigned int**)malloc(MAX_DEGREE * sizeof(vector unsigned int*));
	polyn->coeff_size = (int*) malloc(MAX_DEGREE * sizeof(int));
	polyn->degree = MAX_DEGREE;

	for(i=0; i<MAX_DEGREE; i++)
	{
		polyn->coeff[i] = (vector unsigned int*)malloc(MAX_SIZE * sizeof(vector unsigned int));
	}

	initVector(polyn);
}

inline void initPoly2(poly* polyn, int degree)
{
	int i __attribute__ ((aligned(16)));

	polyn->coeff = (vector unsigned int**)malloc(degree * sizeof(vector unsigned int*));
	polyn->coeff_size = (int*) malloc(degree * sizeof(int));
	polyn->degree = degree;

	for(i=0; i<degree; i++)
	{
		polyn->coeff[i] = (vector unsigned int*)malloc(MAX_SIZE * sizeof(vector unsigned int));
	}

	initVector(polyn);
}

/**
 * The entry point of the program
 * @param  	spe_id the id of the SPU
 * @param	parm pointer to the control structure in the ea
 * @return	0 for no error
 */
int main(unsigned long long spe_id, unsigned long long parm)
{
	int tag_id = 0;
	(void)spe_id;
	(void)parm;

	unsigned int myInput;
	unsigned int myRint;
	unsigned int from;
	unsigned int til;
	unsigned int i;
	unsigned int isPrime;
	int kk=0;

	/*transfer the control structure from the ea into ls*/
	spu_writech(MFC_WrTagMask, -1);
	mfc_get(&cb, (unsigned long) parm, sizeof(control_block), tag_id, 0, 0);
	mfc_read_tag_status_all();

	/*initialize the variables needed*/
	myRint = cb.r;
	from = cb.from;
	til = cb.until;

	/*build the vector that is containing the number to test*/
	vector unsigned int myInputV = (vector unsigned int){cb.mpm_num[0], cb.mpm_num[1], cb.mpm_num[2], cb.mpm_num[3]};

	/*the following lines are computing the log2 of the number*/
	vector unsigned int sumv = spu_cntlz(myInputV);
	myLimit = 0;
	/*to be added properly*/
	for(kk=0; kk<4; kk++)
	{
		int x = spu_extract(sumv, kk);
		myLimit += x;
		if(x != 32)
		{
			break;
		}
	}
	myLimit = 127 - myLimit;

	basisPoly = (poly*) malloc(sizeof(poly));
	poly* resultPoly = (poly*) malloc(sizeof(poly));;

	/*initialize the polynoms*/
	initPoly2(basisPoly, 2);
	initPoly2(resultPoly, myRint);

	basisPoly->degree = 2;
	basisPoly->coeff[1][0] = (vector unsigned int){0, 0, 0, 1};
	basisPoly->coeff[0][0] = (vector unsigned int){0, 0, 0, (unsigned int)1};

	rightPoly_degree = (modR) + 1;
	rightPoly_coeffmodr = (vector unsigned int){0, 0, 0, 1};
	rightPoly_coeff0 = (vector unsigned int){0, 0, 0, (unsigned int)1};


	dummyPoly = (poly*)malloc(sizeof(poly));
	initPoly2(dummyPoly, myRint);

	/*the modulo for the coefficients*/
	vector unsigned int rV = (vector unsigned int){0, 0, 0, myRint};
	mpm_mod(&rV, &myInputV, 1, &rV, 1);
	modR = spu_extract(rV, 3);

	/*flag for the result...1=prime*/
	isPrime = 1;

	/*the exponentiation number*/
	expv = myInputV;
	expv_size = 1;

	/*allocate memory for two temporal MPM numbers used by the calculations*/
	temporal = (vector unsigned int*)malloc(2 * MAX_SIZE * sizeof(vector unsigned int));
	temporal_size = 1;
	temporal2 = (vector unsigned int*)malloc(2 * MAX_SIZE * sizeof(vector unsigned int));
	temporal2_size = 1;

	for(i=from; i<= til; i++) {
		basisPoly->coeff[0][0] = (vector unsigned int){0, 0, 0, (unsigned int)i};
		rightPoly_coeff0 = (vector unsigned int){0, 0, 0, (unsigned int)i};

		polyExp_Mod_int(resultPoly, myRint);
		/*
		 * check if (x + a)^n mod (x^r - 1) != x^(n mod r) + a
		 */

		/*if the polynoms are not equal means that the number being tested is not prime*/
		if(__builtin_expect(polyIsEqual(resultPoly, &rightPoly_coeffmodr, &rightPoly_coeff0) != 1, 1)) {
			/*set the result flag to false*/
			isPrime=0;
			/*inform the PPU about the sad news*/
			spu_write_out_mbox(isPrime);
			break;
		}
	}

	free(temporal);
	free(temporal2);

	/*inform the PPU that the program has finished*/
	spu_write_out_mbox(EXIT_MESSAGE);

	return 0;
}

/**
 * Prints out an instance of the poly-struct
 * @param  	the poly that should be printed
 * @return      void
 */
void polyPrint(poly* in) {
	int i __attribute__ ((aligned(16)));
	(void) printf("%d (", in->degree);
	for (i=in->degree-1; i>=0; i--) {
		if(mpm_cmpeq2( in->coeff[i],  in->coeff_size[i], &incNr, incNr_size) == TRUE_MPM){
			(void) printf(" +x^%d ", i);
		} else {
			if(mpm_cmpeq2( in->coeff[i],  in->coeff_size[i], &zeroNr, zeroNr_size) == FALSE_MPM) {
(void) printf(" +[%u,%u,%u,%u]x^%d ", spu_extract(*in->coeff[i],0), spu_extract(*in->coeff[i],1), spu_extract(*in->coeff[i],2), spu_extract(*in->coeff[i],3), i);
			}
		}

	}
	(void) printf(")\n");
}

/*
 * Makes a deep copy of the input struct and return it
 *
 * @param  	the polynoms that should be multiplied
 * @return      the result polynom
 */
void polyCopy(poly* resultPoly, poly* inPoly) {
	int i __attribute__ ((aligned(16)));
	int j __attribute__ ((aligned(16)));

	initVector(resultPoly);

	resultPoly->degree = inPoly->degree;
	for (i=0; i<= (inPoly->degree-1); i++) {
		resultPoly->coeff_size[i] = inPoly->coeff_size[i];
		/*for(j=0; j<MAX_SIZE; j++)
		{
			resultPoly->coeff[i][j] = inPoly->coeff[i][j];
		}*/
		resultPoly->coeff[i][0] = inPoly->coeff[i][0];
	}
}

/**
 * Tests the equality of a result polynom and the basis polynom
 * @param  	the poly that should be printed
 * @param	the higher degree coefficient of the basis polynom
 # @param	the lower degree coefficient of the basis polynom
 * @return      1 if the polynoms are equal, 0 otherwise
 */
int polyIsEqual(poly* left, vector unsigned int* right_coeff1, vector unsigned int* right_coeff0) {
	int i __attribute__ ((aligned(16)));

	/*are the needed highest coefficients equal?*/
	if (mpm_cmpeq2(left->coeff[rightPoly_degree - 1], left->coeff_size[rightPoly_degree - 1], right_coeff1, 1) == FALSE_MPM)
	{
		/*no they are not*/
		return 0;
	}

	/*are the needed lowest coefficients equal?*/
	if(mpm_cmpeq2(left->coeff[0], left->coeff_size[0], right_coeff0, 1) == FALSE_MPM)
	{
		/*no they are not*/
		return 0;
	}

	/*test the rest of the polynom to have coefficients equal to zero*/
	for(i=1; __builtin_expect(i<left->degree, 1); i++) {
		/*except the on already tested*/
		if(__builtin_expect(i == rightPoly_degree-1, 0))
		{
			i++;
		}
		if(__builtin_expect((mpm_cmpeq2( left->coeff[i],  left->coeff_size[i], &zeroNr, zeroNr_size) == FALSE_MPM), 1))
		{
			/*different means bad news again*/
			return 0;
		}
	}

	return 1;
}

/**
 * This function is multiplying a polynom with a binom
 * @param  	the poly that should be printed
 * @warning 	do not use the method as it is for this program not properly created
 * @return      void
 */
void polyMultBinom(poly* poly, vector unsigned int *up, vector unsigned int* down)
{
	/*
	Probleme:
		1. wird nicht mod r gemacht fur die position
		2. up wird jetzt nicht betrachtet und wird default 1 genommen
	*/
	int i __attribute__ ((aligned(16)));
	i = 0;
	/*because of the binom degree is only one greater*/
	poly->degree++;

	/*carry the first coefficient*/
	poly->coeff[poly->degree - 1][0] = poly->coeff[poly->degree - 2][0];

	/*for all the rest coefficients*/
	for(i = poly->degree - 2; i > 0; i--)
	{
		/*compute the product with the coefficient low degree*/
		temporal_size = 2;
		mpm_mul(temporal, down, 1, poly->coeff[i], poly->coeff_size[i]);
		temporal_size = mpm_add2(temporal, temporal, temporal_size, poly->coeff[i-1], poly->coeff_size[i-1]);

		/*mod the result*/
		mpm_mod(poly->coeff[i], temporal, temporal_size, &expv, expv_size);
		poly->coeff_size[i] = expv_size;
	}

	/*the last coefficient is without add*/
	temporal_size = 2;
	mpm_mul(temporal, down, 1, poly->coeff[0], poly->coeff_size[0]);
	mpm_mod(poly->coeff[0], temporal, temporal_size, &expv, expv_size);
	poly->coeff_size[0] = expv_size;
}

/**
 * Exp function  for a given polynom
 * @param  	the poly that should be exponentiated
 * @param 	the power
 * @param	the modulo of the coeffients
 * @return      void
 */
void polyExp_Mod_int(poly* resultPoly, int r){
	int k __attribute__ ((aligned(16)));

	/*initialize the result polynom using the basis polynom*/
	initVector(resultPoly);
	resultPoly->degree = basisPoly->degree;
	resultPoly->coeff[basisPoly->degree-1][0] = basisPoly->coeff[1][0];
	resultPoly->coeff[0][0] = basisPoly->coeff[0][0];

	/*repeated squaring and multiplying*/
	for(k = myLimit - 1; k > -1; k--){
		polyMult_Mod(dummyPoly, resultPoly, resultPoly, r);/*quadrat*/

		if(testBitVector(&expv, expv_size, k) == 1){
			polyMult_Mod(resultPoly, dummyPoly, basisPoly, r);/*multiply with basis*/
		} else {
			polyCopy(resultPoly, dummyPoly); /*simple copy*/
		}
	}
}

/**
 * Multiplication of two given polynoms
 * @param  	the poly that should be printed
 * @param	the other polynom
 * @param	the modulo for the coefficients
 * @return      void
 */
void polyMult_Mod(poly* resultPoly, poly* left, poly* right, int r) {
	int i 		__attribute__ ((aligned(16)));
	int j 		__attribute__ ((aligned(16)));
	int position	__attribute__ ((aligned(16)));

	initVector(resultPoly);

	/*compute the new degree of the resulting polynom*/
	if((left->degree + right->degree - 1) < r){
		resultPoly->degree = left->degree + right->degree - 1;
	} else {
		resultPoly->degree = r;
	}

	/*compute the coefficients*/
	for (i=0; __builtin_expect(i<= (left->degree-1), 1); i++) {
		for (j=0; __builtin_expect(j<= (right->degree-1), 1); j++) {
			position = (i+j)%r;

			/*multiply the coefficients*/
			temporal_size = left->coeff_size[i] + right->coeff_size[j];
			mpm_mul(temporal, left->coeff[i], left->coeff_size[i], right->coeff[j], right->coeff_size[j]);
			/*added the result to the existing coefficient*/
			temporal_size = mpm_add2(temporal, temporal, temporal_size, resultPoly->coeff[position], resultPoly->coeff_size[position]);

			/*mod the result and store it in the resulting polynom*/
			mpm_div(temporal2, resultPoly->coeff[position], temporal, temporal_size, &expv, expv_size);
			resultPoly->coeff_size[position] = expv_size;
		}
	}
}

/**
 * Get the bit of a given number at a given position
 * @param	the number to operate on
 * @param	size of the number
 * @param	the position to get the bit
 * return     	void
 */
int testBitVector(vector unsigned int* number, int number_size, int bit){
	/*function is using intrinsics and only suitable for 128 bit numbers*/
	/*promote the position to a vector*/
	vector unsigned int divv = spu_promote((unsigned int)bit, 3);
	/*right shift the number*/
	vector unsigned int divvc = spu_rlmask(divv, -5);
	/*shift left*/
	vector unsigned int divvr = spu_sub(divv, spu_sl(divvc, 5));

	/*the 0 component of the array because this version is only for 128 bit numbers*/
	vector unsigned int shifted = spu_rlmaskqwbyte(number[0], -spu_extract(divvc, 3));
	/*right shift*/
	shifted = spu_rlmask(shifted, -spu_extract(divvr, 3));
	/*and the results*/
	shifted = spu_and(shifted, incNr);

	/*extract the bit*/
	return spu_extract(shifted, 3);
}
