
#include "libcoin.h"

SEXP R_LinstatExpCov (const SEXP data, const SEXP inputs, 
                      const SEXP y, const SEXP weights) {

    SEXP ans, p;
    int i, *iinputs, *thisweights;
       
    thisweights = Calloc(LENGTH(weights), int);
    iinputs = LOGICAL(inputs);
    PROTECT(ans = allocVector(VECSXP, LENGTH(data)));
    for (i = 0; i < LENGTH(data); i++) {
        if (iinputs[i]) {
            SET_VECTOR_ELT(ans, i, p = allocVector(VECSXP, 4));
            C_LinstatExpCov(VECTOR_ELT(data, i), y, weights, thisweights, p);
        }
    }
    UNPROTECT(1);
    Free(thisweights);
    return(ans);
}
