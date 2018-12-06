/*
 * @file newton_gmp.c
 * This is a root-finding algorithm which uses the newton-
 * iteration
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */


#include <stdbool.h>
#include "gmp.h"

/**
 * Returns true if myInput has a perfect power or false if not. 
 * for further information about the newton-iteration see 
 * http://mathworld.wolfram.com/NewtonsMethod.html or
 * http://mathworld.wolfram.com/NewtonsIteration.html
 * 
 * @param	myInput		an integer value
 * @return      		true or false
 */
int newton_it(mpz_t a)
{
	int i, j, myResult;

	mpz_t xn;
	mpz_t myOldxn;
	mpz_t dummy1, dummy2;

	i=2;

	myResult = false;
	
	/* initialising the mpz-var. */
	mpz_init(myOldxn);
	mpz_init(dummy1);
	mpz_init(dummy2);
	mpz_init(xn);

	mpz_set_ui(dummy1, 2);
	mpz_cdiv_q(xn, a, dummy1);

	while (mpz_cmp_ui(xn,2)>0)
	{
		mpz_set(myOldxn, xn);

		for (j = 1; j < 1000; j++)
		{


			/* those next steps implements the same as
			 * xn = xn - (pow(xn, i) - myInput) / (i * pow(xn, i-1));
			 * but using the gmp-library
			 */
			mpz_pow_ui(dummy1, xn, i);
			mpz_pow_ui(dummy2, xn, i-1);
			mpz_sub(dummy1, dummy1, a);
			mpz_mul_ui(dummy2, dummy2, i);
			mpz_cdiv_q(dummy1, dummy1, dummy2);
			mpz_sub(xn, xn, dummy1);

			/* if the new calculated xn is bigger or equal the last one,
			 * then break
			 */
			if (mpz_cmp(xn,myOldxn)>=0)
			{
				mpz_set(xn, myOldxn);
				break;
			}
			mpz_set(myOldxn, xn);
		}
		mpz_pow_ui(dummy1, xn, i);

		/* if xn^i == a then a is a perfect power !!! */
		if (mpz_cmp(dummy1, a)==0)
		{
			myResult = true;
		}
		i++;

	}

	return myResult;
}
