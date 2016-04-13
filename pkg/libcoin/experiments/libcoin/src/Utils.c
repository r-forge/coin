
#include "libcoin.h"
#include "helpers.h"

/**
    C- and R-interface to La_svd 
    *\param jobu
    *\param jobv
    *\param x
    *\param s
    *\param u
    *\param v
    *\param method
*/

void CR_La_svd(int dim, char *jobu, char *jobv, double *x, double *s, double *u, double *v) {

    int *xdims, n, p, lwork, info = 0;
    int *iwork;
    double *work, *xvals, tmp;

    xvals = Calloc(dim * dim, double);
    /* work on a copy of x */
    Memcpy(xvals, x, (size_t) (dim * dim));

    iwork= (int *) Calloc(8*dim, int);

    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgesdd)(jobu,
                     &dim, &dim, xvals, &dim, s,
		     u, &dim,
		     v, &dim,
		     &tmp, &lwork, iwork, &info);
    if (info != 0)
        error(("error code %d from Lapack routine '%s'"), info, "dgesdd");
    lwork = (int) tmp;
    work = Calloc(lwork, double);
    F77_CALL(dgesdd)(jobu, 
                     &dim, &dim, xvals, &dim, s,
                     u, &dim,
		     v, &dim,
		     work, &lwork, iwork, &info);
    if (info != 0)
        error(("error code %d from Lapack routine '%s'"), info, "dgesdd");
    Free(work); Free(xvals); Free(iwork);
}

/**
    C-interface to CR_La_svd
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

void CR_svd (SEXP x, SEXP s, SEXP u, SEXP v) {

    int P;
    double *du, *dv, *ds;

    if (!isMatrix(x) || !isReal(x))
        error("x is not a real matrix");

    du = REAL(u);
    dv = REAL(v);
    ds = REAL(s);
    P = NCOL(x);
    if (NROW(x) != P) error("svd p x error");

    for (int p = 0; p < P*P; p++) {
        if (p < P) ds[p] = 0.0;
        du[p] = 0.0;
        dv[p] = 0.0;
    }
    
    CR_La_svd(P, "S", "", REAL(x), ds, du, dv);
}


/**
    R-interface to CR_La_svd
    *\param x matrix
    *\param svdmem an object of class `svd_mem'
*/

SEXP R_svd (SEXP x) {

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
    Moore-Penrose inverse of a matrix
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
    *\param ans return value; an object of class `ExpectCovarMPinv'
*/

void C_MPinv (SEXP x, double tol, SEXP s, SEXP u, SEXP v, double *MPinv, double *rank) {

    int i, j, P, k, *positive;
    double *ds, *du, *dvt;
    
    CR_svd(x, s, u, v);
    ds = REAL(s);
    du = REAL(u);
    dvt = REAL(v);
    /* this may be the reduced dimension! Use the first p elements only!!!*/
    P = NCOL(x);

    if (tol * ds[0] > tol) tol = tol * ds[0];

    positive = Calloc(P, int); 
    
    rank[0] = 0.0;
    for (i = 0; i < P; i++) {
        if (ds[i] > tol) {
            positive[i] = 1;
            rank[0] += 1.0;
        } 
    }
    
    for (j = 0; j < P; j++) {
        if (positive[j]) {
            for (i = 0; i < P; i++)
                du[j * P + i] *= (1 / ds[j]);
        }
    }
    
    for (i = 0; i < P; i++) {
        for (j = 0; j < P; j++) {
            MPinv[j * P + i] = 0.0;
            for (k = 0; k < P; k++) {
                if (positive[k])
                    MPinv[j * P + i] += dvt[i * P + k] * du[P * k + j];
            }
        }
    }

    Free(positive);
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
    SET_VECTOR_ELT(ans, 1, rank = allocVector(REALSXP, 1));

    C_MPinv(x, REAL(tol)[0], s, u, v, REAL(MP), REAL(rank));
    
    UNPROTECT(1);
    return(ans);
}
