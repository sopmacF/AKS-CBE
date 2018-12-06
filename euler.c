/*
 * @file euler.c
 * the euler-function (PI(n)), it returns the number of 
 * positive integers less than or equal to n that are coprime to n
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "euler.h"

/*
 #define _DEBUG 
 */

/**
 * Returns the PI(n) using a division-test to prime factor n
 * for further information about the PI-Function see 
 * http://mathworld.wolfram.com/TotientFunction.html or
 * http://en.wikipedia.org/wiki/Euler's_totient_function
 * 
 * The factoring method is based on a trial division
 * 
 * @param	myInput		an integer value
 * @return      		PI(myInput) coprime to myInput
 */
int euler_probDiv(int myInput) {

	float myResult = 1;
	int mySaveInput = myInput;
	int myPower = 0; /* info/debug-var. */
	int myDivisor= false;
	int myCurDivisor= false; /* info/debug-var. */

	int myMax=(int) (myInput/2);

	/* try to divide myInput by 2 as often as possible (without remainder)*/
	while (myInput%2==0) {
		myInput=myInput/2;
		myPower+=1;
		if (!myCurDivisor) {
			/* calculates the distinct prime factors once per factor */
			myResult = 0.5; /* the same as result = result * (1.0 - (1.0/2)); */
			myDivisor = true;
			myCurDivisor = true;
		}
	}
#ifdef _DEBUG
	if(myCurDivisor) {
		(void) printf("2 ^ %d\n", myPower);
		(void) printf("result = %.2f\n", myResult);
	}
#endif

	myPower=0;

	myCurDivisor = false;

	/* try to divide myInput by 3 as often as possible (without remainder)*/
	while (myInput%3==0) {
		myInput=myInput/3;
		myPower+=1;
		if (!myCurDivisor) {
			/* calculates the distinct prime factors once per factor */
			myResult = myResult * (1.0 - (1.0/3));
			myDivisor = true;
			myCurDivisor = true;
		}
	}
#ifdef _DEBUG
	if(myCurDivisor) {
		(void) printf("3 ^ %d\n", myPower);
		(void) printf("result = %.2f\n", myResult);
	}
#endif

	myPower=0;

	int t=5;
	int d=2;

	myMax = (myInput<myMax) ? myInput : myMax;

	myCurDivisor = false;

	/* 
	 * try to divide myInput by t as often as possible; t  => (n(6) +/- 1)	
	 */
	while (t<=myMax) {
		if ((t>sqrt(mySaveInput)) && (!myDivisor)) {
#ifdef _DEBUG
			(void) printf("%d seems to be prime\n", mySaveInput);
#endif
			break;
		}
		while (myInput%t==0) {
			myInput=myInput/t;
			myPower+=1;
			if (!myCurDivisor) {
				/* calculates the distinct prime factors once per factor */
				myResult = myResult * (1.0 - (1.0/t));
				myDivisor = true;
				myCurDivisor = true;
			}
		}
#ifdef _DEBUG
		if(myCurDivisor) {
			(void) printf("%d ^ %d\n", t, myPower);
			(void) printf("result = %.2f\n", myResult);
		}
#endif
		t=t+d;
		d=6-d;

		myCurDivisor = false;
		myPower=0;
	}

	if (!myDivisor) {
		/* myInput seems to be a prime number so PI(myInput) = myInput -1 */
		myResult = mySaveInput - 1;
	} else {
		myResult = myResult * mySaveInput;
	}

	return myResult;
}
