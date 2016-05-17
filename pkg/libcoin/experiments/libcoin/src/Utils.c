
#include "libcoin_internal.h"
#include "Tables.h"
#include "Sums.h"

int NLEVELS(SEXP x) 
{
   return(LENGTH(getAttrib(x, R_LevelsSymbol)));
}

int NROW(SEXP x) 
{
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(LENGTH(x));
    return(INTEGER(a)[0]);
}
    
int NCOL(SEXP x) 
{
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    return(INTEGER(a)[1]);
}

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker(const double *A, const int m, const int n,
                 const double *B, const int r, const int s, int overwrite,
                 double *ans)
{
    int i, j, k, l, mr, js, ir;
    double y;

    if (overwrite) {
        for (i = 0; i < m * r * n * s; i++) ans[i] = 0.0;
    }

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++)
                    ans[(js + l) * mr + ir + k] += y * B[l * r + k];
            }
        }
    }
}  

void C_kronecker_sym(const double *A, const int m, 
                     const double *B, const int r, int overwrite,
                     double *ans)
{
    int i, j, k, l, mr, js, ir, s, n, tmp, mrns;
    double y;

    mr = m * r;
    n = m;
    s = r;

    if (overwrite) {
        for (i = 0; i < mr * (mr + 1) / 2; i++) ans[i] = 0.0;
    }

    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j <= i; j++) {
            js = j * s;
            y = A[S(i, j, m)]; 
            for (k = 0; k < r; k++) {
                for (l = 0; l < (j < i ? s : k + 1); l++) {
                    ans[S(ir + k, js + l, mr)] += y * B[S(k, l, r)]; 
                }
            }
        }
    }
}  


/* MP inv of symmetric matrix in lower triangular packed form */

void C_MPinv_sym (double *x, int n, double tol, double *dMP, int *rank) {

    SEXP ans;
    double *val, *vec, dtol, *rx, *work, valinv;
    int valzero = 0, info = 0, kn;

    rx = Calloc(n * (n + 1) / 2, double);
    Memcpy(rx, x, n * (n + 1) / 2);
    work = Calloc(3 * n, double);
    val = Calloc(n, double);
    vec = Calloc(n * n, double);
    
    F77_CALL(dspev)("V", "L", &n, rx, val, vec, &n, work,
                    &info);
                                            
    dtol = val[n - 1] * tol;

    for (int k = 0; k < n; k++)
        valzero += (val[k] < dtol); 
    rank[0] = n - valzero;

    for (int i = 0; i < n * (n + 1) / 2; i++) dMP[i] = 0.0;
    
    for (int k = valzero; k < n; k++) {
        valinv = 1 / val[k];
        kn = k * n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                /* MP is symmetric */
                dMP[S(i, j, n)] += valinv * vec[kn + i] * vec[kn + j];
            }
        }
    }
    Free(rx); Free(work); Free(val); Free(vec);
}
