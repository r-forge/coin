
#include "libcoin.h"
#include "Utils.h"
#include "helpers.h"

/**
    R-interface to CR_La_svd
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_svd (SEXP x) 
{
    int P;
    SEXP ans, s, u, v;
    
    P = NROW(x);
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, s = allocVector(REALSXP, P));
    SET_VECTOR_ELT(ans, 1, u = allocMatrix(REALSXP, P, P));
    SET_VECTOR_ELT(ans, 2, v = allocMatrix(REALSXP, P, P));

    CR_svd(x, s, u, v);
    
    UNPROTECT(1);
    return(ans);
}


/**
    R-interface to C_MPinv 
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_MPinv (SEXP x, SEXP tol) {

    SEXP ans, MP, rank, s, u, v;
    int P;

    if (!isMatrix(x) || !isReal(x))
        error("R_MPinv: x is not a real matrix");

    if (NROW(x) != NCOL(x)) 
        error("R_MPinv: x is not a square matrix");

    if (!isReal(tol) || LENGTH(tol) != 1)
        error("R_MPinv: tol is not a scalar real");

    P = NROW(x);
    /* potentially, the effective dimension was reduced
    if (p != INTEGER(GET_SLOT(svdmem, PL2_pSym))[0])
        error("R_MPinv: dimensions don't match");
    */

    s = allocVector(REALSXP, P);
    u = allocMatrix(REALSXP, P, P);
    v = allocMatrix(REALSXP, P, P);

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, MP = allocMatrix(REALSXP, P, P));
    SET_VECTOR_ELT(ans, 1, rank = allocVector(INTSXP, 1));

    C_MPinv(x, REAL(tol)[0], s, u, v, REAL(MP), INTEGER(rank));
    
    UNPROTECT(1);
    return(ans);
}
