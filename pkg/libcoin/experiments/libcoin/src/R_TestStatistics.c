
#include "libcoin_internal.h"
#include "TestStatistics.h"
#include "R_LinearStatistic.h"
#include "Utils.h"
#include "Contrasts.h"
#include "MemoryAccess.h"

SEXP R_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = CR_quadform(LinearStatistic, Expectation, MPinv);
    UNPROTECT(1);
    return(ans);
}

SEXP R_maxtype(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP type, SEXP tol)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    if (INTEGER(type)[0] == ALTERNATIVE_twosided)
        REAL(ans)[0] = CR_maxabsstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_less)
        REAL(ans)[0] = CR_maxstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_greater)
        REAL(ans)[0] = CR_minstat(LinearStatistic, Expectation, CoVariance, tol);
    UNPROTECT(1);
    return(ans);
}

SEXP R_orderedmaxsel(SEXP LECV, SEXP teststat, SEXP tol, SEXP minbucket)
{
    SEXP ans, wmax, maxstat;
    double *contrasts, *V, *ExpX, xtab, total;
    int P, Q, mb, nc = 0, start = 0, stop = 0;

    P = C_get_P(LECV);
    Q = C_get_Q(LECV);
    mb = INTEGER(minbucket)[0];
    ExpX = C_get_ExpectationX(LECV);

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, wmax = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, maxstat = allocVector(REALSXP, 1));

    xtab = 0.0;
    for (int p = 0; p < P; p++) {
        xtab += ExpX[p];
        if (xtab > mb) {
            start = p;
            break;
        }
    }
    xtab = 0.0;
    for (int p = (P - 1); p >= 0; p--) {
        xtab += ExpX[p];
        if (xtab > mb) {
            stop = p;
            break;
        }
    }
    if (start >= stop)
        error("cannot find admissible split");
    nc = stop - start;

    contrasts = Calloc(P * nc, double);
    
    for (int p = 0; p < P * nc; p++) contrasts[p] = 0.0;

    for (int p = start; p < stop; p++) {
        for (int pp = 0; pp <= p; pp++)
            contrasts[pp + (p - start) * P] = 1.0;
    }

    if (INTEGER(teststat)[0] == 1) {                            
        if (C_get_varonly(LECV)) {
            V = C_get_Variance(LECV);
        } else {
            V = C_get_Covariance(LECV);
        }
        C_contrasts_marginal_maxabsstat(C_get_LinearStatistic(LECV),
                                        C_get_Expectation(LECV), 
                                        V,
                                        contrasts, P, Q, 
                                        nc,
                                        REAL(tol)[0], 
                                        INTEGER(wmax), REAL(maxstat));
    } else {
        if (C_get_varonly(LECV))
            error("cannot compute quadratic form from variance only");
        C_contrasts_marginal_quadform(C_get_LinearStatistic(LECV),
                                      C_get_Expectation(LECV), 
                                      C_get_Covariance(LECV),
                                      contrasts, P, Q, 
                                      nc,
                                      REAL(tol)[0], 
                                      INTEGER(wmax), REAL(maxstat));
    }

    INTEGER(wmax)[0] += start;

    Free(contrasts);
    UNPROTECT(1);
    return(ans);
}                                      

SEXP R_unorderedmaxsel(SEXP LECV, SEXP teststat, SEXP tol, SEXP minbucket)
{
    SEXP ans, wmax, maxstat, index;
    double *contrasts, *V, *ExpX, xtab, total, *indl;
    int P, Q, mb, nc = 0, start = 0, stop = 0;

    P = C_get_P(LECV);
    if (P >= 31)
        error("cannot search for unordered splits in >= 31 levels");
    Q = C_get_Q(LECV);
    mb = INTEGER(minbucket)[0];
    ExpX = C_get_ExpectationX(LECV);

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, wmax = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, maxstat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, P));

    total = 0.0;
    for (int p = 0; p < P; p++) total += ExpX[p];
              
    /* number of possible binary splits */
    int mi = 1;  
    for (int l = 1; l < P; l++) mi *= 2;
    contrasts = Calloc(mi * P, double);

    nc = 0;
    for (int j = 1; j < mi; j++) { /* go though all splits */
        indl = contrasts + P * nc;
         
         /* indl determines if level p is left or right */
         int jj = j;
         for (int l = 1; l < P; l++) {
             indl[l] = (jj%2);
             jj /= 2;
         }
                
         xtab = 0.0;
         for (int p = 0; p < P; p++)
             xtab += indl[p] * ExpX[p];
         if (xtab > mb && (total - xtab) > mb)
             nc++;
    }

    if (INTEGER(teststat)[0] == 1) {                            
        if (C_get_varonly(LECV)) {
            V = C_get_Variance(LECV);
        } else {
            V = C_get_Covariance(LECV);
        }
        C_contrasts_marginal_maxabsstat(C_get_LinearStatistic(LECV),
                                        C_get_Expectation(LECV), 
                                        V,
                                        contrasts, P, Q, 
                                        nc,
                                        REAL(tol)[0], 
                                        INTEGER(wmax), REAL(maxstat));
    } else {
        if (C_get_varonly(LECV))
            error("cannot compute quadratic form from variance only");
        C_contrasts_marginal_quadform(C_get_LinearStatistic(LECV),
                                      C_get_Expectation(LECV), 
                                      C_get_Covariance(LECV),
                                      contrasts, P, Q, 
                                      nc,
                                      REAL(tol)[0], 
                                      INTEGER(wmax), REAL(maxstat));
    }
    
    for (int p = 0; p < P; p++)
        INTEGER(index)[p] = (int) contrasts[p + INTEGER(wmax)[0] * P] + 1;
    
    Free(contrasts);
    UNPROTECT(1);
    return(ans);
}                                      

