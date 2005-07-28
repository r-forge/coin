
/**
    S4 classes 
    *\file Classes.c
    *\author $Author$
    *\date $Date$
*/

#include "CI_common.h"

SEXP 
    CI_expectationSym,
    CI_covarianceSym,
    CI_sumweightsSym;

SEXP coin_init(void) {
    CI_expectationSym = install("expectation");
    CI_covarianceSym = install("covariance");
    CI_sumweightsSym = install("sumweights");
    return(R_NilValue);
}
