
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

void CR_La_svd(int dim, char *jobu, char *jobv, double *x, 
               double *s, double *u, double *v) 
{
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

void CR_svd (SEXP x, SEXP s, SEXP u, SEXP v) 
{
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
    Moore-Penrose inverse of a matrix
    *\param x matrix
    *\param tol a tolerance bound
    *\param svdmem an object of class `svd_mem'
    *\param ans return value; an object of class `ExpectCovarMPinv'
*/

void C_MPinv (SEXP x, double tol, SEXP s, SEXP u, SEXP v, double *MPinv, 
              int *rank)
{
    int i, j, P, k;
    double *ds, *du, *dvt;
    
    CR_svd(x, s, u, v);
    ds = REAL(s);
    du = REAL(u);
    dvt = REAL(v);
    P = NCOL(x);

    if (tol * ds[0] > tol) tol = tol * ds[0];

    rank[0] = 0;
    for (i = 0; i < P; i++) {
        if (ds[i] > tol) rank[0]++;
    }
    
    for (j = 0; j < P; j++) {
        if (ds[j] > tol) {
            for (i = 0; i < P; i++)
                du[j * P + i] *= (1 / ds[j]);
        }
    }
    
    for (i = 0; i < P; i++) {
        for (j = 0; j < P; j++) {
            MPinv[j * P + i] = 0.0;
            for (k = 0; k < P; k++) {
                if (ds[k] > tol)
                    MPinv[j * P + i] += dvt[i * P + k] * du[P * k + j];
            }
        }
    }
}

/* MP inv of symmetric matrix in lower triangular packed form */

SEXP R_MPinv_sym (SEXP x, SEXP tol) {

    SEXP ans, rank, MP;
    double *val, *vec, *dMP, dtol, *rx, *work, valinv;
    int n, valzero = 0, info = 0, kn;

    n = (int) (.5 + sqrt(.25 + 2 * LENGTH(x))) - 1;
    
    rx = Calloc(LENGTH(x), double);
    Memcpy(rx, REAL(x), LENGTH(x));
    work = Calloc(3 * n, double);
    val = Calloc(n, double);
    vec = Calloc(n * n, double);
    
    F77_CALL(dspev)("V", "L", &n, rx, val, vec, &n, work,
                    &info);
                                            
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, MP = allocMatrix(REALSXP, n, n));
    SET_VECTOR_ELT(ans, 1, rank = allocVector(INTSXP, 1));
    dMP = REAL(MP);
    
    dtol = val[n - 1] * REAL(tol)[0];

    for (int k = 0; k < n; k++)
        valzero += (val[k] < dtol); 
    INTEGER(rank)[0] = n - valzero;

    for (int i = 0; i < n * n; i++) dMP[i] = 0.0;
    
    for (int k = valzero; k < n; k++) {
        valinv = 1 / val[k];
        kn = k * n;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                /* MP is symmetric */
                dMP[j * n + i] += valinv * vec[kn + i] * vec[kn + j];
                dMP[i * n + j] = dMP[j * n + i];
            }
        }
    }
    Free(rx); Free(work); Free(val); Free(vec);
    UNPROTECT(1);
    return(ans);
}
