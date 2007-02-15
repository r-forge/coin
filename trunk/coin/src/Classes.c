
/**
    S4 classes 
    *\file Classes.c
    *\author $Author$
    *\date $Date$
*/

#include "coin_common.h"

SEXP 
    coin_expectationSym,
    coin_covarianceSym,
    coin_sumweightsSym;

SEXP coin_init(void) {
    coin_expectationSym = install("expectation");
    coin_covarianceSym = install("covariance");
    coin_sumweightsSym = install("sumweights");
    return(R_NilValue);
}
