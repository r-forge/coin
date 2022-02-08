
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/stats_package.h> /* for S_rcont2 */
#include <R_ext/stats_stubs.h>  /* for S_rcont2 */


static void
C_Smirnov_sim(int nrow, int ncol, const int nrowt[], const int ncolt[], int n,
              int B, int *observed, int twosided, double *fact, int *jwork, double *results)
{
    /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
    fact[0] = fact[1] = 0.;
    for(int i = 2; i <= n; i++)
        fact[i] = fact[i - 1] + log(i);

    GetRNGstate();

    for(int iter = 0; iter < B; ++iter) {
        S_rcont2(nrow, ncol, nrowt, ncolt, n, fact, jwork, observed);
        double S = 0., diff = 0.;
        int cs0 = 0, cs1 = 0;
        for (int j = 0; j < nrow; j++) {
            cs0 += observed[j];
            cs1 += observed[nrow + j];
            diff = ((double) cs0) / ncolt[0] - ((double) cs1) / ncolt[1];
            if (twosided) diff = fabs(diff);
            if (diff > S) S = diff;
        }
        results[iter] = S;
    }

    PutRNGstate();

    return;
}

SEXP Smirnov_sim(SEXP sr, SEXP sc, SEXP sB, SEXP twosided)
{
    sr = PROTECT(coerceVector(sr, INTSXP));
    sc = PROTECT(coerceVector(sc, INTSXP));
    int nr = LENGTH(sr), nc = LENGTH(sc), B = asInteger(sB);
    if (nc != 2)
        error("Smirnov statistic only defined for two groups"); 
    int n = 0, *isr = INTEGER(sr);
    for (int i = 0; i < nr; i++) {
        /* avoid integer overflow */
        if (n > INT_MAX - isr[i]) 
            error("Sample size too large");
        n += isr[i];
    }
    int *observed = (int *) R_alloc(nr * nc, sizeof(int));
    double *fact = (double *) R_alloc(n+1, sizeof(double));
    int *jwork = (int *) R_alloc(nc, sizeof(int));
    SEXP ans = PROTECT(allocVector(REALSXP, B));
    C_Smirnov_sim(nr, nc, isr, INTEGER(sc), n, B, observed, 
                  INTEGER(twosided)[0], fact, jwork, REAL(ans));
    UNPROTECT(3);
    return ans;
}
