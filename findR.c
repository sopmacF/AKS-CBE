/*
 * @file findR.c
 * this function returns the smallest r
 * such that such that or() > log^2 n
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */

#include <stdbool.h>
#include "findR.h"
#include "gmp.h"

/**
 * Returns the smallest r
 * such that such that or() > log^2 n
 * 
 * @param	n			the input n
 * @return      		the r-value found
 */
int findR(mpz_t n) {
	mpz_t r;
	mpz_t myLimit;
	mpz_t myLogN;
	mpz_t myResult;

	bool failed = false;
	int value = 0;

	mpz_init(r);
	mpz_init(myLimit);
	mpz_init(myResult);
	mpz_init(myLogN);

	mpz_set_ui(r, 2);

	mpz_set_ui(myLogN, (int)mpz_sizeinbase(n, 2));

	mpz_mul(myLimit, myLogN, myLogN);

	//mpz_mul_ui(myLimit, myLimit, 4);

	while ((mpz_cmp(n, r)>0)) {
		if (mpz_probab_prime_p(r, 5)) {
			mpz_t myCounter;

			mpz_init(myCounter);

			mpz_set_ui(myCounter, 1);
			failed=false;

			while ((mpz_cmp(myLimit, myCounter)>=0) && (!failed)) {
				mpz_set_ui(myResult, 0);

				mpz_powm(myResult, n, myCounter, r);

				if (mpz_cmp_ui(myResult, 1)==0) {

					failed = true;
				}
				mpz_add_ui(myCounter, myCounter, 1);
			}
			if (!failed) {
				/* the smallest r found !!! */
				return mpz_get_ui(r);
			}
		}
		mpz_add_ui(r, r, 1);
	}
	if (mpz_cmp(n, r)==0) {
		/* composite !!! */
		value = -1;

	}
	return value;
}
