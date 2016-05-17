
#include "libcoin_internal.h"
#include "Utils.h"

SEXP R_MPinv_sym (SEXP x, SEXP tol) {

    SEXP ans, MP, rank;
    int n;
    
    n = (int) (.5 + sqrt(.25 + 2 * LENGTH(x))) - 1;

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, MP = allocVector(REALSXP, n * (n + 1) / 2));
    SET_VECTOR_ELT(ans, 1, rank = allocVector(INTSXP, 1));
            
    C_MPinv_sym(REAL(x), n, REAL(tol)[0], REAL(MP), INTEGER(rank));
    
    UNPROTECT(1);
    return(ans);
}

    