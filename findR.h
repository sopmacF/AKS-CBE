/*
 * @file findR.h
 * this function returns the smallest r
 * such that such that or() > log^2 n
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */

#include <stdbool.h>
#include "gmp.h"

#ifndef FINDR_H_
#define FINDR_H_

int findR(mpz_t n);

#endif /*FINDR_H_*/
