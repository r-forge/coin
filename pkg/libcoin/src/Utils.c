
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

    if (n == 1) {
        if (x[0] > tol) {
            dMP[0] = 1 / x[0];
            rank[0] = 1;
        } else {
            dMP[0] = 0;
            rank[0] = 0;
        }
    } else {
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
}

void
rcont2(int *nrow, int *ncol,
       int *nrowt, int *ncolt, int *ntotal,
       double *fact, int *jwork, int *matrix)
{
    int j, l, m, ia, ib, ic, jc, id, ie, ii, nll, nlm, nr_1, nc_1;
    double x, y, dummy, sumprb;
    Rboolean lsm, lsp;

    nr_1 = *nrow - 1;
    nc_1 = *ncol - 1;

    ib = 0; /* -Wall */

    /* Construct random matrix */
    for (j = 0; j < nc_1; ++j)
	jwork[j] = ncolt[j];

    jc = *ntotal;

    for (l = 0; l < nr_1; ++l) { /* -----  matrix[ l, * ] ----- */
	ia = nrowt[l];
	ic = jc;
	jc -= ia;/* = n_tot - sum(nr[0:l]) */

	for (m = 0; m < nc_1; ++m) {
	    id = jwork[m];
	    ie = ic;
	    ic -= id;
	    ib = ie - ia;
	    ii = ib - id;

	    if (ie == 0) { /* Row [l,] is full, fill rest with zero entries */
		for (j = m; j < nc_1; ++j)
		    matrix[l + j * *nrow] = 0;
		ia = 0;
		break;
	    }

	    /* Generate pseudo-random number */
	    dummy = unif_rand();

	    do {/* Outer Loop */

		/* Compute conditional expected value of MATRIX(L, M) */

		nlm = (int)(ia * (id / (double) ie) + 0.5);
		x = exp(fact[ia] + fact[ib] + fact[ic] + fact[id]
			- fact[ie] - fact[nlm]
			- fact[id - nlm] - fact[ia - nlm] - fact[ii + nlm]);
		if (x >= dummy)
		    break;
		if (x == 0.)/* MM: I haven't seen this anymore */
		    error("rcont2 [%d,%d]: exp underflow to 0; algorithm failure", l, m);

		sumprb = x;
		y = x;
		nll = nlm;

		do {
		    /* Increment entry in row L, column M */
		    j = (int)((id - nlm) * (double)(ia - nlm));
		    lsp = (j == 0);
		    if (!lsp) {
			++nlm;
			x = x * j / ((double) nlm * (ii + nlm));
			sumprb += x;
			if (sumprb >= dummy)
			    goto L160;
		    }

		    do {
			R_CheckUserInterrupt();

			/* Decrement entry in row L, column M */
			j = (int)(nll * (double)(ii + nll));
			lsm = (j == 0);
			if (!lsm) {
			    --nll;
			    y = y * j / ((double) (id - nll) * (ia - nll));
			    sumprb += y;
			    if (sumprb >= dummy) {
				nlm = nll;
				goto L160;
			    }
			    /* else */
			    if (!lsp)
				break;/* to while (!lsp) */
			}
		    } while (!lsm);

		} while (!lsp);

		dummy = sumprb * unif_rand();

	    } while (1);

L160:
	    matrix[l + m * *nrow] = nlm;
	    ia -= nlm;
	    jwork[m] -= nlm;
	}
	matrix[l + nc_1 * *nrow] = ia;/* last column in row l */
    }

    /* Compute entries in last row of MATRIX */
    for (m = 0; m < nc_1; ++m)
	matrix[nr_1 + m * *nrow] = jwork[m];

    matrix[nr_1 + nc_1 * *nrow] = ib - matrix[nr_1 + (nc_1-1) * *nrow];

    return;
}
