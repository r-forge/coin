
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "Distributions.h"
#include "MemoryAccess.h"
#include "Contrasts.h"

SEXP R_ChisqTest(SEXP LEV, SEXP linstat, SEXP tol, SEXP lower, SEXP give_log) 
{
    SEXP ans, stat, pval;
    double *MPinv, *pv, st, *ls, *ex;
    int rank, P, Q, PQ, B;
    
    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;

    if (C_get_varonly(LEV))
        error("cannot compute quadratic form based on variances only");

    MPinv = C_get_MPinv(LEV);
    C_MPinv_sym(C_get_Covariance(LEV), PQ, REAL(tol)[0], MPinv, &rank);
        
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    
    REAL(stat)[0] = C_quadform(PQ, C_get_LinearStatistic(LEV),
                               C_get_Expectation(LEV), MPinv);
    if (LENGTH(linstat) == 0) {
        REAL(pval)[0] = C_chisq_pvalue(REAL(stat)[0], rank, INTEGER(lower)[0],
                                       INTEGER(give_log)[0]);
    } else {
        B = NCOL(linstat);
        ls = REAL(linstat);
        pv = REAL(pval);
        st = REAL(stat)[0];
        ex = C_get_Expectation(LEV);
        pv[0] = 0.0;
        for (int i = 0; i < B; i++) {
            if (C_quadform(PQ, ls + PQ * i, ex, MPinv)> st)
                pv[0] = pv[0] + 1.0;
        }
        if (INTEGER(give_log)[0]) {
            if (INTEGER(lower)[0]) {
                pv[0] = log1p(- pv[0] / B);
            } else {
                pv[0] = log(pv[0]) - log(B);
            }
        } else {
            if (INTEGER(lower)[0]) {
                pv[0] = 1 - pv[0] / B;
            } else {
                pv[0] = pv[0] / B;
            }
        }
    }

    UNPROTECT(1);
    return(ans);
}

SEXP R_MaxabsstatTest(SEXP LEV, SEXP linstat, SEXP tol, SEXP lower, 
                      SEXP give_log, SEXP maxpts, SEXP releps, SEXP abseps)
{
    SEXP ans, stat, pval;
    double st, *ex, *cv, *ls, tl, *pv;
    int P, Q, PQ, B;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
            
    if (C_get_varonly(LEV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));

    REAL(stat)[0] =  C_maxabsstat_Covariance(PQ, C_get_LinearStatistic(LEV), 
                                             C_get_Expectation(LEV), 
                                             C_get_Covariance(LEV), REAL(tol)[0]);
   
    if (LENGTH(linstat) == 0) {
        REAL(pval)[0] = C_maxabsstat_pvalue(REAL(stat)[0], C_get_Covariance(LEV),
                                            PQ,
                                            INTEGER(maxpts)[0], REAL(releps)[0], 
                                            REAL(abseps)[0], REAL(tol)[0]);
        if (INTEGER(give_log)[0]) {
            if (INTEGER(lower)[0]) {
                REAL(pval)[0] = log1p(- REAL(pval)[0]);
            } else {
                REAL(pval)[0] = log(REAL(pval)[0]);
            }
        } else {
            if (INTEGER(lower)[0])
                REAL(pval)[0] = 1 - REAL(pval)[0];
        }
    } else {
        B = NCOL(linstat);
        pv = REAL(pval);
        st = REAL(stat)[0];
        ls = REAL(linstat);
        ex = C_get_Expectation(LEV);
        cv = C_get_Covariance(LEV);
        tl = REAL(tol)[0];
        REAL(pval)[0] = 0.0;
        for (int i = 0; i < B; i++) {
            if (C_maxabsstat_Covariance(PQ, ls + PQ * i, ex, cv, tl) > st)
                pv[0] = pv[0] + 1.0;
        }
        if (INTEGER(give_log)[0]) {
            if (INTEGER(lower)[0]) {
                pv[0] = log1p(- pv[0] / B);
            } else {
                pv[0] = log(pv[0]) - log(B);
            }
        } else {
            if (INTEGER(lower)[0]) {
                pv[0] = 1 - pv[0] / B;
            } else {
                pv[0] = pv[0] / B;
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}

SEXP R_MaxstatTest_ordered(SEXP LEV, SEXP linstat, SEXP teststat, SEXP tol, 
                           SEXP minbucket, SEXP lower, SEXP give_log)
{
    SEXP ans, index, stat, pval;
    double *contrasts, *V, *ExpX, xtab, tmp, *pv, *ls, st;
    int P, Q, PQ, B, mb, nc = 0, start = 0, stop = 0, itmp;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
    mb = INTEGER(minbucket)[0];
    ExpX = C_get_ExpectationX(LEV);

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, 1));

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
        if (C_get_varonly(LEV)) {
            V = C_get_Variance(LEV);
        } else {
            V = C_get_Covariance(LEV);
        }
        C_contrasts_marginal_maxabsstat(C_get_LinearStatistic(LEV),
                                        C_get_Expectation(LEV), 
                                        V,
                                        contrasts, P, Q, 
                                        nc,
                                        REAL(tol)[0], 
                                        INTEGER(index), REAL(stat));
    } else {
        if (C_get_varonly(LEV))
            error("cannot compute quadratic form from variance only");
        C_contrasts_marginal_quadform(C_get_LinearStatistic(LEV),
                                      C_get_Expectation(LEV), 
                                      C_get_Covariance(LEV),
                                      contrasts, P, Q, 
                                      nc,
                                      REAL(tol)[0], 
                                      INTEGER(index), REAL(stat));
    }

    INTEGER(index)[0] += start;

    if (LENGTH(linstat) > 0) {
        st = REAL(stat)[0];
        pv = REAL(pval);
        ls = REAL(linstat);
        B = NCOL(linstat);

        for (int i = 0; i < B; i++) {
            if (INTEGER(teststat)[0] == 1) {                            
                C_contrasts_marginal_maxabsstat(ls + PQ * i,
                                                C_get_Expectation(LEV), 
                                                V,
                                                contrasts, P, Q, 
                                                nc,
                                                REAL(tol)[0], 
                                                &itmp, &tmp);
            } else {
                C_contrasts_marginal_quadform(ls + PQ * i,
                                              C_get_Expectation(LEV), 
                                              C_get_Covariance(LEV),
                                              contrasts, P, Q, 
                                              nc,
                                              REAL(tol)[0], 
                                              &itmp, &tmp);
            }
            if (tmp > st) pv[0] = pv[0] + 1.0;
        }
        if (INTEGER(give_log)[0]) {
            if (INTEGER(lower)[0]) {
                pv[0] = log1p(- pv[0] / B);
            } else {
                pv[0] = log(pv[0]) - log(B);
            }
        } else {
            if (INTEGER(lower)[0]) {
                pv[0] = 1 - pv[0] / B;
            } else {
                pv[0] = pv[0] / B;
            }
        }
    }
    Free(contrasts);
    UNPROTECT(1);
    return(ans);
}                                      

SEXP R_MaxstatTest_unordered(SEXP LEV, SEXP linstat, SEXP teststat, SEXP tol, 
                             SEXP minbucket, SEXP lower, SEXP give_log)
{
    SEXP ans, stat, index, pval;
    double *contrasts, *V, *ExpX, xtab, total, *indl, *ls, *pv, tmp, st;
    int P, Q, PQ, B, mb, nc = 0, wmax;

    P = C_get_P(LEV);
    if (P >= 31)
        error("cannot search for unordered splits in >= 31 levels");
    Q = C_get_Q(LEV);
    PQ = P * Q;
    mb = INTEGER(minbucket)[0];
    ExpX = C_get_ExpectationX(LEV);

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
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
        if (C_get_varonly(LEV)) {
            V = C_get_Variance(LEV);
        } else {
            V = C_get_Covariance(LEV);
        }
        C_contrasts_marginal_maxabsstat(C_get_LinearStatistic(LEV),
                                        C_get_Expectation(LEV), 
                                        V,
                                        contrasts, P, Q, 
                                        nc,
                                        REAL(tol)[0], 
                                        &wmax, REAL(stat));
    } else {
        if (C_get_varonly(LEV))
            error("cannot compute quadratic form from variance only");
        C_contrasts_marginal_quadform(C_get_LinearStatistic(LEV),
                                      C_get_Expectation(LEV), 
                                      C_get_Covariance(LEV),
                                      contrasts, P, Q, 
                                      nc,
                                      REAL(tol)[0], 
                                      &wmax, REAL(stat));
    }
    
    for (int p = 0; p < P; p++)
        INTEGER(index)[p] = (int) contrasts[p + wmax * P] + 1;

    if (LENGTH(linstat) > 0) {
        st = REAL(stat)[0];
        pv = REAL(pval);
        ls = REAL(linstat);
        B = NCOL(linstat);

        for (int i = 0; i < B; i++) {
            if (INTEGER(teststat)[0] == 1) {                            
                C_contrasts_marginal_maxabsstat(ls + PQ * i,
                                                C_get_Expectation(LEV), 
                                                V,
                                                contrasts, P, Q, 
                                                nc,
                                                REAL(tol)[0], 
                                                &wmax, &tmp);
            } else {
            C_contrasts_marginal_quadform(ls + PQ * i,
                                          C_get_Expectation(LEV), 
                                          C_get_Covariance(LEV),
                                          contrasts, P, Q, 
                                          nc,
                                          REAL(tol)[0], 
                                          &wmax, &tmp);
            }
            if (tmp > st) pv[0] = pv[0] + 1.0;
        }
        if (INTEGER(give_log)[0]) {
            if (INTEGER(lower)[0]) {
                pv[0] = log1p(- pv[0] / B);
            } else {
                pv[0] = log(pv[0]) - log(B);
            }
        } else {
            if (INTEGER(lower)[0]) {
                pv[0] = 1 - pv[0] / B;
            } else {
                pv[0] = pv[0] / B;
            }
        }
    }
    Free(contrasts);
    UNPROTECT(1);
    return(ans);
}                                      
