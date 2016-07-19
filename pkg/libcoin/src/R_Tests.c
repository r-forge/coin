
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

SEXP R_MaxtypeTest(SEXP LEV, SEXP linstat, SEXP tol, SEXP alternative, SEXP lower, 
                   SEXP give_log, SEXP maxpts, SEXP releps, SEXP abseps)
{
    SEXP ans, stat, pval;
    double st, *ex, *cv, *ls, tl, *pv;
    int P, Q, PQ, B, vo, alt;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
            
    if (C_get_varonly(LEV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));

    REAL(stat)[0] =  C_maxtype(PQ, C_get_LinearStatistic(LEV), 
                                   C_get_Expectation(LEV), 
                                   C_get_Covariance(LEV), 
                                   C_get_varonly(LEV),
                                   REAL(tol)[0],
                                   INTEGER(alternative)[0]);

    if (LENGTH(linstat) == 0) {
        if (C_get_varonly(LEV) && PQ > 1) {
            REAL(pval)[0] = NA_REAL;
            UNPROTECT(1);
            return(ans);
        }
        REAL(pval)[0] = C_maxtype_pvalue(REAL(stat)[0], C_get_Covariance(LEV),
                                         PQ, INTEGER(alternative)[0], INTEGER(lower)[0],
                                         INTEGER(give_log)[0],
                                         INTEGER(maxpts)[0], REAL(releps)[0], 
                                         REAL(abseps)[0], REAL(tol)[0]);
    } else {
        B = NCOL(linstat);
        pv = REAL(pval);
        st = REAL(stat)[0];
        ls = REAL(linstat);
        ex = C_get_Expectation(LEV);
        cv = C_get_Covariance(LEV);
        vo = C_get_varonly(LEV);
        alt = INTEGER(alternative)[0];
        tl = REAL(tol)[0];
        REAL(pval)[0] = 0.0;
        for (int i = 0; i < B; i++) {
            if (alt == ALTERNATIVE_less) {
                if (LE(C_maxtype(PQ, ls + PQ * i, ex, cv, vo, tl, alt), st, tl))
                    pv[0] = pv[0] + 1.0;
            } else {
                if (GE(C_maxtype(PQ, ls + PQ * i, ex, cv, vo, tl, alt), st, tl))
                    pv[0] = pv[0] + 1.0;
            }
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
    double tmp, *pv, *ls, st;
    int P, Q, PQ, B, mb, itmp;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
    mb = INTEGER(minbucket)[0];

    if (C_get_varonly(LEV))
        error("cannot maximally selected statistics from variance only");

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, 1));

    if (INTEGER(teststat)[0] == 1) {                            
        C_ordered_maxabsstand_Xfactor(C_get_LinearStatistic(LEV),
                               C_get_Expectation(LEV), 
                               C_get_Covariance(LEV),
                               P, Q, 
                               C_get_ExpectationX(LEV),
                               mb,
                               REAL(tol)[0], 
                               INTEGER(index), REAL(stat));
    } else {
        C_ordered_quadform_Xfactor(C_get_LinearStatistic(LEV),
                             C_get_Expectation(LEV),
                             C_get_Covariance(LEV),
                             P, Q,
                             C_get_ExpectationX(LEV),
                             mb,
                             REAL(tol)[0],
                             INTEGER(index), REAL(stat));
    }

    /* no admissible split found */
    if (INTEGER(index)[0] < 0) {
        REAL(stat)[0] = 0.0;
        REAL(pval)[0] = 1.0;
        INTEGER(index)[0] = NA_INTEGER;
        UNPROTECT(1);
        return(ans);
    } else {
        INTEGER(index)[0]++; /* R indexing */
    }

    if (LENGTH(linstat) > 0) {
        st = REAL(stat)[0];
        pv = REAL(pval);
        ls = REAL(linstat);
        B = NCOL(linstat);

        for (int i = 0; i < B; i++) {
            if (INTEGER(teststat)[0] == 1) {                            
                C_ordered_maxabsstand_Xfactor(ls + PQ * i,
                                       C_get_Expectation(LEV),
                                       C_get_Covariance(LEV),
                                       P, Q,
                                       C_get_ExpectationX(LEV),
                                       mb,
                                       REAL(tol)[0],
                                       &itmp, &tmp);
            } else {
                C_ordered_quadform_Xfactor(ls + PQ * i,
                                     C_get_Expectation(LEV),
                                     C_get_Covariance(LEV),
                                     P, Q,
                                     C_get_ExpectationX(LEV),
                                     mb,
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
    UNPROTECT(1);
    return(ans);
}                                      

SEXP R_MaxstatTest_unordered(SEXP LEV, SEXP linstat, SEXP teststat, SEXP tol, 
                             SEXP minbucket, SEXP lower, SEXP give_log)
{
    SEXP ans, stat, index, pval;
    double *contrasts, *ExpX, sumleft, totalsum, *indl, *ls, *pv, tmp, st;
    int P, Pnonzero, Q, PQ, B, mb, nc = 0, wmax, *levels;;


    if (C_get_varonly(LEV))
        error("cannot compute maximally selected statistics from variance only");

    P = C_get_P(LEV);
    ExpX = C_get_ExpectationX(LEV);
    sumleft = 0.0;
    totalsum = 0.0;
    Pnonzero = 0;
    for (int p = 0; p < P; p++)  {
        totalsum += ExpX[p];
        if (ExpX[p] > 0) Pnonzero++;
    }

    levels = Calloc(Pnonzero, int);
    nc = 0;
    for (int p = 0; p < P; p++) {
        if (ExpX[p] > 0) {
            levels[nc] = p;
            nc++;
        }
    }
    
    if (Pnonzero >= 31)
        error("cannot search for unordered splits in >= 31 levels");
    Q = C_get_Q(LEV);
    PQ = P * Q;
    mb = INTEGER(minbucket)[0];

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, index = allocVector(INTSXP, P));

    for (int p = 0; p < P; p++) {
        if (ExpX[p] == 0.0) INTEGER(index)[p] = NA_INTEGER;
    }
              
    /* number of possible binary splits */
    int mi = 1;  
    for (int l = 1; l < Pnonzero; l++) mi *= 2;
    contrasts = Calloc(mi * P, double);
    for (int j = 0; j < mi * P; j++) contrasts[j] = 0.0;

    indl = Calloc(Pnonzero, double);
    for (int p = 0; p < Pnonzero; p++) indl[p] = 0.0;

    nc = 0;
    for (int j = 1; j < mi; j++) { /* go though all splits */
         
         /* indl determines if level p is left or right */
         int jj = j;
         for (int l = 1; l < Pnonzero; l++) {
             indl[l] = (jj%2);
             jj /= 2;
         }
                
         sumleft = 0.0;
         for (int p = 0; p < Pnonzero; p++)
             sumleft += indl[p] * ExpX[levels[p]];

         if (sumleft > mb && (totalsum - sumleft) > mb) {
             for (int p = 0; p < Pnonzero; p++)
                 contrasts[P * nc + levels[p]] = indl[p];
             nc++;
         }
    }

    if (INTEGER(teststat)[0] == 1) {                            
        C_contrasts_marginal_maxabsstand(C_get_LinearStatistic(LEV),
                                        C_get_Expectation(LEV), 
                                        C_get_Covariance(LEV),
                                        contrasts, P, Q, 
                                        nc,
                                        REAL(tol)[0], 
                                        &wmax, REAL(stat));
    } else {
        C_contrasts_marginal_quadform(C_get_LinearStatistic(LEV),
                                      C_get_Expectation(LEV), 
                                      C_get_Covariance(LEV),
                                      contrasts, P, Q, 
                                      nc,
                                      REAL(tol)[0], 
                                      &wmax, REAL(stat));
    }
    
    /* no admissible split found */
    if (wmax < 0) {
        REAL(stat)[0] = 0.0;
        REAL(pval)[0] = 1.0;
        /* INDEX == NA anyhow */
        UNPROTECT(1);
        return(ans);
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
                C_contrasts_marginal_maxabsstand(ls + PQ * i,
                                                C_get_Expectation(LEV), 
                                                C_get_Covariance(LEV),
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
    Free(contrasts); Free(indl); Free(levels);
    UNPROTECT(1);
    return(ans);
}                                      
