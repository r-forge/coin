
#include "libcoin_internal.h"
#include "Tables.h"
#include "Sums.h"

int C_nlevels (SEXP x) 
{
   return(LENGTH(getAttrib(x, R_LevelsSymbol)));
}

int NROW (SEXP x) 
{
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(LENGTH(x));
    return(INTEGER(a)[0]);
}
    
int NCOL (SEXP x) 
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

void C_kronecker (const double *A, const int m, const int n,
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

void C_kronecker_sym (const double *A, const int m, 
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
