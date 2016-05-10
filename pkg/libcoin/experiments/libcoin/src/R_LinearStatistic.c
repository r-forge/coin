
#include "libcoin.h"
#include "LinearStatistic.h"
#include "helpers.h"
#include "Sums.h"
#include "Tables.h"
#include "Distributions.h"

SEXP R_LinearStatistic(SEXP x, SEXP y)
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    C_LinearStatistic(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_weights(SEXP x, SEXP y, SEXP weights)
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    C_LinearStatistic_weights(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), 
                              INTEGER(weights), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_subset(SEXP x, SEXP y, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    C_LinearStatistic_subset(REAL(x), NROW(x), NCOL(x), REAL(y), NCOL(y), 
                             INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_weights_subset(SEXP x, SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    C_LinearStatistic_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                     REAL(y), NCOL(y), INTEGER(weights), 
                                     INTEGER(subset), LENGTH(subset), 
                                     REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_LinearStatistic_2d(SEXP x, SEXP y, SEXP table) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(REALSXP, NCOL(x), NCOL(y)));
    C_LinearStatistic_2d(REAL(x), NROW(x), NCOL(x), REAL(y), NROW(y), NCOL(y), 
                         INTEGER(table), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence(SEXP y) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence(REAL(y), NROW(y), NCOL(y), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_weights(SEXP y, SEXP weights) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_weights(REAL(y), NROW(y), NCOL(y), 
                                   INTEGER(weights), 
                                   C_sum(INTEGER(weights), LENGTH(weights)), 
                                   REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_subset(SEXP y, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_subset(REAL(y), NROW(y), NCOL(y), 
                                  INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationInfluence_weights_subset(SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_ExpectationInfluence_weights_subset(REAL(y), NROW(y), NCOL(y), 
                                          INTEGER(weights), 
                                          C_sum_subset(INTEGER(weights), 
                                                       LENGTH(weights), 
                                                       INTEGER(subset), 
                                                       LENGTH(subset)),
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceInfluence(SEXP y) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence(y));
    PROTECT(ans = allocVector(REALSXP, NCOL(y) * (NCOL(y) + 1) / 2));
    C_CovarianceInfluence(REAL(y), NROW(y), NCOL(y), REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_weights(SEXP y, SEXP weights) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights(y, weights));
    PROTECT(ans = allocVector(REALSXP, NCOL(y) * (NCOL(y) + 1) / 2));
    C_CovarianceInfluence_weights(REAL(y), NROW(y), NCOL(y), 
                                  INTEGER(weights), 
                                  C_sum(INTEGER(weights), LENGTH(weights)), 
                                  REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_subset(SEXP y, SEXP subset) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_subset(y, subset));
    PROTECT(ans = allocVector(REALSXP, NCOL(y) * (NCOL(y) + 1) / 2));
    C_CovarianceInfluence_subset(REAL(y), NROW(y), NCOL(y), 
                                 INTEGER(subset), LENGTH(subset), 
                                 REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_CovarianceInfluence_weights_subset(SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights_subset(y, weights, subset));
    PROTECT(ans = allocVector(REALSXP, NCOL(y) * (NCOL(y) + 1) / 2));
    C_CovarianceInfluence_weights_subset(REAL(y), NROW(y), NCOL(y), 
                                         INTEGER(weights), 
                                         C_sum_subset(INTEGER(weights), 
                                                      LENGTH(weights), 
                                                      INTEGER(subset), 
                                                      LENGTH(subset)),
                                         INTEGER(subset), LENGTH(subset), 
                                         REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_VarianceInfluence(SEXP y) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence(y));
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_VarianceInfluence(REAL(y), NROW(y), NCOL(y), REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_VarianceInfluence_weights(SEXP y, SEXP weights) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights(y, weights));
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_VarianceInfluence_weights(REAL(y), NROW(y), NCOL(y), 
                                INTEGER(weights), 
                                C_sum(INTEGER(weights), LENGTH(weights)), 
                                REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_VarianceInfluence_subset(SEXP y, SEXP subset) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_subset(y, subset));
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_VarianceInfluence_subset(REAL(y), NROW(y), NCOL(y), 
                               INTEGER(subset), LENGTH(subset), 
                               REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_VarianceInfluence_weights_subset(SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans, ExpInf;
    
    PROTECT(ExpInf = R_ExpectationInfluence_weights_subset(y, weights, subset));
    PROTECT(ans = allocVector(REALSXP, NCOL(y)));
    C_VarianceInfluence_weights_subset(REAL(y), NROW(y), NCOL(y), 
                                       INTEGER(weights), 
                                       C_sum_subset(INTEGER(weights), 
                                                    LENGTH(weights), 
                                                    INTEGER(subset), 
                                                    LENGTH(subset)),
                                       INTEGER(subset), LENGTH(subset), 
                                       REAL(ExpInf), REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_ExpectationX(SEXP x) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX(REAL(x), NROW(x), NCOL(x), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_weights(SEXP x, SEXP weights) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_weights(REAL(x), NROW(x), NCOL(x), 
                           INTEGER(weights), 
                           REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_subset(SEXP x, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_subset(REAL(x), NROW(x), NCOL(x), 
                          INTEGER(subset), LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationX_weights_subset(SEXP x, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_ExpectationX_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                  INTEGER(weights), 
                                  INTEGER(subset), LENGTH(subset), 
                                  REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX(SEXP x) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x) * (NCOL(x) + 1) / 2));
    C_CovarianceX(REAL(x), NROW(x), NCOL(x), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_weights(SEXP x, SEXP weights) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x) * (NCOL(x) + 1) / 2));
    C_CovarianceX_weights(REAL(x), NROW(x), NCOL(x), 
                          INTEGER(weights), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_subset(SEXP x, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x) * (NCOL(x) + 1) / 2));
    C_CovarianceX_subset(REAL(x), NROW(x), NCOL(x), 
                                 INTEGER(subset), LENGTH(subset), 
                                 REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceX_weights_subset(SEXP x, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x) * (NCOL(x) + 1) / 2));
    C_CovarianceX_weights_subset(REAL(x), NROW(x), NCOL(x), 
                                          INTEGER(weights), 
                                          INTEGER(subset), LENGTH(subset), 
                                          REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_VarianceX(SEXP x) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_VarianceX(REAL(x), NROW(x), NCOL(x), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_VarianceX_weights(SEXP x, SEXP weights) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_VarianceX_weights(REAL(x), NROW(x), NCOL(x), 
                                   INTEGER(weights), 
                                   REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_VarianceX_subset(SEXP x, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_VarianceX_subset(REAL(x), NROW(x), NCOL(x), 
                       INTEGER(subset), LENGTH(subset), 
                       REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_VarianceX_weights_subset(SEXP x, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, NCOL(x)));
    C_VarianceX_weights_subset(REAL(x), NROW(x), NCOL(x), 
                               INTEGER(weights), INTEGER(subset), 
                              LENGTH(subset), REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_ExpectationLinearStatistic(SEXP ExpInf, SEXP ExpX) 
{
    SEXP ans;
    int P, Q;
    
    P = LENGTH(ExpX);
    Q = LENGTH(ExpInf);

    PROTECT(ans = allocMatrix(REALSXP, P, Q));
    C_ExpectationLinearStatistic(P, Q, REAL(ExpInf), REAL(ExpX), 0, REAL(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_CovarianceLinearStatistic(SEXP CovInf, SEXP ExpX, SEXP CovX, 
                                 SEXP sumweights) 
{
    SEXP ans, PP_tmp;
    
    int P = LENGTH(ExpX);
    int Q = (int) (.5 + sqrt(.25 + 2 * LENGTH(CovInf))) - 1;
    int PQ = P * Q;
    PROTECT(ans = allocVector(REALSXP, PQ * (PQ + 1) / 2));
    PROTECT(PP_tmp = allocVector(REALSXP, P * (P + 1) / 2));
    C_CovarianceLinearStatistic(P, Q, REAL(CovInf), REAL(ExpX), REAL(CovX), 
                                INTEGER(sumweights)[0], REAL(PP_tmp), 0, REAL(ans));
    UNPROTECT(2);
    return(ans);
}

SEXP R_VarianceLinearStatistic(SEXP VarInf, SEXP ExpX, SEXP VarX, SEXP sumweights) 
{
    SEXP ans, P_tmp;
    
    int P = LENGTH(ExpX);
    int Q = LENGTH(VarInf);
    PROTECT(ans = allocVector(REALSXP, P * Q));
    PROTECT(P_tmp = allocVector(REALSXP, P));
    C_VarianceLinearStatistic(P, Q, REAL(VarInf), REAL(ExpX), REAL(VarX), 
                              INTEGER(sumweights)[0], REAL(P_tmp), 0, REAL(ans));
    UNPROTECT(2);
    return(ans);
}


SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset, SEXP block, SEXP varonly)

{
    SEXP ans, L, E, V; 
    int N, P, Q, Nlevel, *sumweights, *table, *subset_tmp, tmp;
    double *ExpInf, *CovInf, *work;
    
    N = NROW(x);
    P = NCOL(x);
    Q = NCOL(y);
    
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, L = allocVector(REALSXP, P * Q));
    SET_VECTOR_ELT(ans, 1, E = allocVector(REALSXP, P * Q));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, P * Q));
    } else  {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, P * Q * (P * Q + 1) / 2));
    }

    if (LENGTH(block) == 0) {
        table = Calloc(1, int);
        table[0] = LENGTH(subset);
        sumweights = Calloc(1, int);
        if (LENGTH(weights) == 0) {
            sumweights[0] = 0;
        } else {
            if (LENGTH(subset) == 0) {
                sumweights[0] = C_sum(INTEGER(weights), LENGTH(weights));
            } else {
                sumweights[0] = C_sum_subset(INTEGER(weights), LENGTH(weights), 
                                             INTEGER(subset), LENGTH(subset));
            }
        }
        ExpInf = Calloc(Q, double);
        if (INTEGER(varonly)[0]) {
           CovInf = Calloc(Q, double);
           work = Calloc(3 * P, double);
        } else {
            CovInf = Calloc(Q * (Q + 1) / 2, double);
            work = Calloc(P + 2 * P * (P + 1) / 2, double);
        }
        C_LinearStatistic_(REAL(x), N, P, REAL(y), Q, INTEGER(weights), sumweights, 
                          INTEGER(subset), table, 1, REAL(L));
        if (INTEGER(varonly)[0]) {
            C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
                sumweights, INTEGER(subset), table, 1, 1, ExpInf, CovInf);
            C_ExpectationVarianceLinearStatistic(REAL(x), N, P, Q, INTEGER(weights),
                sumweights, INTEGER(subset), table, 1, ExpInf, CovInf, work, REAL(E), REAL(V)); 
        } else {
            C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
                sumweights, INTEGER(subset), table, 1, 0, ExpInf, CovInf);
            C_ExpectationCovarianceLinearStatistic(REAL(x), N, P, Q, INTEGER(weights),
                sumweights, INTEGER(subset), table, 1, ExpInf, CovInf, work, REAL(E), REAL(V)); 
        }
        Free(table); Free(sumweights); Free(table); Free(ExpInf); Free(CovInf); Free(work);
    } else {
        Nlevel = C_nlevels(block);
        table = Calloc(Nlevel + 1, int);

        if (LENGTH(subset) == 0) {
            C_1dtable(INTEGER(block), Nlevel + 1, N, table);
            subset_tmp = Calloc(N, int);
            C_setup_subset(N, subset_tmp);
            C_order_wrt_block(subset_tmp, N, INTEGER(block), table, Nlevel + 1);
        } else {
            C_1dtable_subset(INTEGER(block), Nlevel + 1, INTEGER(subset), LENGTH(subset), table);
            subset_tmp = Calloc(LENGTH(subset), int);
            Memcpy(subset_tmp, INTEGER(subset), LENGTH(subset));
            C_order_wrt_block(subset_tmp, LENGTH(subset), INTEGER(block), table, Nlevel + 1);
        }

        sumweights = Calloc(Nlevel, int);
        if (LENGTH(weights) == 0) {        
            for (int b = 0; b < Nlevel; b++) sumweights[b] = 0;
        } else {
            tmp = 0;
            for (int b = 0; b < Nlevel; b++) {
                sumweights[b] = C_sum_subset(INTEGER(weights), LENGTH(weights), subset_tmp + tmp, table[b + 1]);
                tmp = tmp + table[b + 1];
            }
        }
        ExpInf = Calloc(Nlevel * Q, double);
        if (INTEGER(varonly)[0]) {
            CovInf = Calloc(Nlevel * Q, double);
            work = Calloc(3 * P, double);
        } else {
            CovInf = Calloc(Nlevel * Q * (Q + 1) / 2, double);
            work = Calloc(P + 2 * P * (P + 1) / 2, double);
        }
        C_LinearStatistic_(REAL(x), N, P, REAL(y), Q, INTEGER(weights), sumweights,
                           subset_tmp, table + 1, Nlevel, REAL(L));
        if (INTEGER(varonly)[0]) {
            C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
                sumweights, subset_tmp, table + 1, Nlevel, 1, ExpInf, CovInf);
            C_ExpectationVarianceLinearStatistic(REAL(x), N, P, Q, INTEGER(weights),
                sumweights, subset_tmp, table + 1, Nlevel, ExpInf, CovInf, work, REAL(E), REAL(V)); 
        } else {
            C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
                sumweights, subset_tmp, table + 1, Nlevel, 0, ExpInf, CovInf);
            C_ExpectationCovarianceLinearStatistic(REAL(x), N, P, Q, INTEGER(weights),
                sumweights, subset_tmp, table + 1, Nlevel, ExpInf, CovInf, work, REAL(E), REAL(V)); 
        }
        Free(subset_tmp); Free(table); Free(sumweights); Free(table); Free(ExpInf); Free(CovInf);
        Free(work);
    }
    UNPROTECT(1);
    return(ans);
}        

SEXP R_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy, SEXP weights, SEXP subset, SEXP block, SEXP varonly)

{
    SEXP ans, L, E, V; 
    int N, P, Q, Nlevel, *table, *btab, *csum, *rsum, *table2d, sw;
    double *ExpInf, *CovInf, *ExpX, *CovX, *work;
    
    N = LENGTH(ix);
    P = NCOL(x);
    Q = NCOL(y);
    Nlevel = 1;
    if (LENGTH(block) > 0)
        Nlevel = C_nlevels(block);

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, L = allocVector(REALSXP, P * Q));
    SET_VECTOR_ELT(ans, 1, E = allocVector(REALSXP, P * Q));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, P * Q));
    } else  {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, P * Q * (P * Q + 1) / 2));
    }

    table = Calloc((C_nlevels(ix) + 1) * (C_nlevels(iy) + 1) * Nlevel, int);
    table2d = Calloc((C_nlevels(ix) + 1) * (C_nlevels(iy) + 1), int);
    csum = Calloc((C_nlevels(iy) + 1), int);
    rsum = Calloc((C_nlevels(ix) + 1), int);
    
    C_2dtable_(ix, iy, weights, subset, block, table);
    
    for (int b = 0; b < Nlevel; b++) {
        if (b == 0) {
            for (int i = 0; i < (C_nlevels(ix) + 1) * (C_nlevels(iy) + 1); i++)
                table2d[i] = 0;
        }
        for (int i = 0; i < (C_nlevels(ix) + 1); i++) {
            for (int j = 0; j < (C_nlevels(iy) + 1); j++)
                table2d[j * (C_nlevels(ix) + 1) + i] += 
                    table[b * (C_nlevels(ix) + 1) * (C_nlevels(iy) + 1) + 
                          j * (C_nlevels(ix) + 1) + i];
        }
    }

    C_LinearStatistic_2d(REAL(x), NROW(x), P, REAL(y), NROW(y), Q, table2d, REAL(L));

    ExpInf = Calloc(Q, double);
    ExpX = Calloc(P, double);
    if (INTEGER(varonly)[0]) {
        CovInf = Calloc(Q, double);
        CovX = Calloc(P, double);
        work = Calloc(P, double);
    } else {
        CovInf = Calloc(Q * (Q + 1) / 2, double);
        CovX = Calloc(P * (P + 1) / 2, double);
        work = Calloc(P * (P + 1) / 2, double);
    }

    for (int b = 0; b < Nlevel; b++) {
        btab = table + (C_nlevels(ix) + 1) * (C_nlevels(iy) + 1) * b;
        C_colSums_i(btab, (C_nlevels(ix) + 1), (C_nlevels(iy) + 1), csum); 
        csum[0] = 0; /* NA */   
        C_rowSums_i(btab, (C_nlevels(ix) + 1), (C_nlevels(iy) + 1), rsum);    
        rsum[0] = 0; /* NA */
        sw = 0;
        for (int i = 1; i < (C_nlevels(ix) + 1); i++) sw += rsum[i];

        C_ExpectationInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf);
        C_ExpectationX_weights(REAL(x), NROW(x), P, rsum, ExpX);
        C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, REAL(E));
        if (INTEGER(varonly)[0]) {
            C_VarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf, CovInf);
            C_VarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            C_VarianceLinearStatistic(P, Q, CovInf, ExpX, CovX, sw, work, b, REAL(V));
        } else {
            C_CovarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf, CovInf);
            C_CovarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            C_CovarianceLinearStatistic(P, Q, CovInf, ExpX, CovX, sw, work, b, REAL(V));
        }
    }
    Free(table); 
    Free(table2d); 
    Free(csum); 
    Free(rsum); 
    Free(ExpInf); 
    Free(ExpX);
    Free(CovInf); 
    Free(CovX); 
    Free(work);
    UNPROTECT(1);
    return(ans);
}        

