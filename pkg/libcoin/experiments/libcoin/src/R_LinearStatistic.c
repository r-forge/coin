
#include "libcoin_internal.h"
#include "LinearStatistic.h"
#include "Utils.h"
#include "Sums.h"
#include "Tables.h"
#include "Distributions.h"

void RC_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, 
                                       SEXP subset, SEXP block, SEXP varonly, 
                                       double *L, double *E, double *V)

{
    int N, P, Q, Nlevel, *sumweights, *table, *subset_tmp, tmp;
    double *ExpInf, *CovInf, *work;
    
    N = NROW(x);
    if (isInteger(x)) {
        P = NLEVELS(x);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);
    
    Nlevel = 1;
    if (LENGTH(block) > 0)
        Nlevel = NLEVELS(block);

    ExpInf = Calloc(Nlevel * Q, double);
    if (INTEGER(varonly)[0]) {
       CovInf = Calloc(Nlevel * Q, double);
       work = Calloc(3 * P + 1, double);
    } else {
        CovInf = Calloc(Nlevel * Q * (Q + 1) / 2, double);
        work = Calloc(P + 2 * P * (P + 1) / 2 + 1, double);
    }
    table = Calloc(Nlevel + 1, int);
    sumweights = Calloc(Nlevel, int);

    if (Nlevel == 1) {
        table[0] = 0;
        table[1] = LENGTH(subset);
        if (LENGTH(weights) == 0) {
            sumweights[0] = 0;
        } else {
            if (LENGTH(subset) == 0) {
                sumweights[0] = C_sum_(INTEGER(weights), LENGTH(weights));
            } else {
                sumweights[0] = C_sum_subset(INTEGER(weights), LENGTH(weights), 
                                             INTEGER(subset), LENGTH(subset));
            }
        }
        subset_tmp = INTEGER(subset);
    } else {
        if (LENGTH(subset) == 0) {
            C_1dtable_(INTEGER(block), Nlevel + 1, N, table);
            subset_tmp = Calloc(N, int);
            C_setup_subset(N, subset_tmp);
            C_order_wrt_block(subset_tmp, N, INTEGER(block), table, Nlevel + 1);
        } else {
            C_1dtable_subset(INTEGER(block), Nlevel + 1, INTEGER(subset), 
                             LENGTH(subset), table);
            subset_tmp = Calloc(LENGTH(subset), int);
            Memcpy(subset_tmp, INTEGER(subset), LENGTH(subset));
            C_order_wrt_block(subset_tmp, LENGTH(subset), INTEGER(block), 
                              table, Nlevel + 1);
        }

        if (LENGTH(weights) == 0) {        
            for (int b = 0; b < Nlevel; b++) sumweights[b] = 0;
        } else {
            tmp = 0;
            for (int b = 0; b < Nlevel; b++) {
                sumweights[b] = C_sum_subset(INTEGER(weights), LENGTH(weights), 
                                             subset_tmp + tmp, table[b + 1]);
                tmp = tmp + table[b + 1];
            }
        }
    }

    C_LinearStatistic(x, N, P, REAL(y), Q, INTEGER(weights), 
                      sumweights, subset_tmp, table + 1, Nlevel, L);

    if (INTEGER(varonly)[0]) {
        C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Nlevel, 1, ExpInf, CovInf);
        C_ExpectationVarianceLinearStatistic(x, N, P, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Nlevel, ExpInf, CovInf, 
            work, E, V); 
    } else {
        C_ExpectationCovarianceInfluence(REAL(y), N, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Nlevel, 0, ExpInf, CovInf);
        C_ExpectationCovarianceLinearStatistic(x, N, P, Q, INTEGER(weights),
            sumweights, subset_tmp, table + 1, Nlevel, ExpInf, CovInf, work, 
            E, V); 
    }
    if (Nlevel > 1) Free(subset_tmp); 
    
    Free(table); Free(sumweights); Free(ExpInf); 
    Free(CovInf); Free(work);
}        


SEXP R_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, SEXP subset, 
                                      SEXP block, SEXP varonly)

{
    SEXP ans, L, E, V; 
    int P, Q, PQ;
    
    if (isInteger(x)) {
        P = NLEVELS(x);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);
    PQ = P * Q;
    
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, L = allocVector(REALSXP, PQ));
    SET_VECTOR_ELT(ans, 1, E = allocVector(REALSXP, PQ));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, PQ));
    } else  {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, PQ * (PQ + 1) / 2));
    }

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, varonly,
                                      REAL(L), REAL(E), REAL(V));
    UNPROTECT(1);
    return(ans);
}



void RC_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy,
                                          SEXP weights, SEXP subset, SEXP block,
                                          SEXP varonly, double *L, double *E, 
                                          double *V)
{

    int N, P, Q, Nlevel, *table, *btab, *csum, *rsum, *table2d, sw, *iix;
    double *ExpInf, *CovInf, *ExpX, *CovX, *work;

    N = LENGTH(ix);
    if (LENGTH(x) == 0) {
        P = NLEVELS(ix);
    } else {
        P = NCOL(x);
    }
    Q = NCOL(y);

    Nlevel = 1;
    if (LENGTH(block) > 0)
        Nlevel = NLEVELS(block);

    table = Calloc((NLEVELS(ix) + 1) * (NLEVELS(iy) + 1) * Nlevel, int);
    table2d = Calloc((NLEVELS(ix) + 1) * (NLEVELS(iy) + 1), int);
    csum = Calloc((NLEVELS(iy) + 1), int);
    rsum = Calloc((NLEVELS(ix) + 1), int);
    
    RC_2dtable(ix, iy, weights, subset, block, table);
    
    for (int i = 0; i < (NLEVELS(ix) + 1) * (NLEVELS(iy) + 1); i++)
        table2d[i] = 0;
    for (int b = 0; b < Nlevel; b++) {
        for (int i = 0; i < (NLEVELS(ix) + 1); i++) {
            for (int j = 0; j < (NLEVELS(iy) + 1); j++)
                table2d[j * (NLEVELS(ix) + 1) + i] += 
                    table[b * (NLEVELS(ix) + 1) * (NLEVELS(iy) + 1) + 
                          j * (NLEVELS(ix) + 1) + i];
        }
    }

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

    if (LENGTH(x) == 0) {
        C_LinearStatistic_2d(ix, LENGTH(ix), P, REAL(y), NROW(y), Q, 
                             table2d, L);
    } else {
        C_LinearStatistic_2d(x, NROW(x), P, REAL(y), NROW(y), Q, 
                             table2d, L);
    }

    for (int b = 0; b < Nlevel; b++) {
        btab = table + (NLEVELS(ix) + 1) * (NLEVELS(iy) + 1) * b;
        C_colSums_i(btab, (NLEVELS(ix) + 1), (NLEVELS(iy) + 1), csum); 
        csum[0] = 0; /* NA */   
        C_rowSums_i(btab, (NLEVELS(ix) + 1), (NLEVELS(iy) + 1), rsum);    
        rsum[0] = 0; /* NA */
        sw = 0;
        for (int i = 1; i < (NLEVELS(ix) + 1); i++) sw += rsum[i];
        C_ExpectationInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf);
        if (LENGTH(x) == 0) {
            for (int p = 0; p < P; p++)
                ExpX[p] = (double) rsum[p + 1];
        } else {
            C_ExpectationX_weights(REAL(x), NROW(x), P, rsum, ExpX);
        }
        C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, E);
        if (INTEGER(varonly)[0]) {
            C_VarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf, 
                                        CovInf);
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P; p++) CovX[p] = ExpX[p];
            } else {
                C_VarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            }
            C_VarianceLinearStatistic(P, Q, CovInf, ExpX, CovX, sw, work, b, V);
        } else {
            C_CovarianceInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf, 
                                          CovInf);
            if (LENGTH(x) == 0) {
                for (int p = 0; p < P * (P + 1) / 2; p++) CovX[p] = 0.0;
                for (int p = 0; p < P; p++) CovX[S(p, p, P)] = ExpX[p];
            } else {
                C_CovarianceX_weights(REAL(x), NROW(x), P, rsum, CovX);
            }
            C_CovarianceLinearStatistic(P, Q, CovInf, ExpX, CovX, sw, work, b, V);
        }
    }
    Free(table); Free(table2d); Free(csum); Free(rsum); Free(ExpInf); Free(ExpX);
    Free(CovInf); Free(CovX); Free(work);
}        


SEXP R_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy, 
                                         SEXP weights, SEXP subset, SEXP block, 
                                         SEXP varonly)
{
    SEXP ans, L, E, V; 
    int P, Q, PQ;
    
    P = NCOL(x);
    Q = NCOL(y);
    PQ = P * Q;

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, L = allocVector(REALSXP, PQ));
    SET_VECTOR_ELT(ans, 1, E = allocVector(REALSXP, PQ));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, PQ));
    } else  {
        SET_VECTOR_ELT(ans, 2, V = allocVector(REALSXP, PQ * (PQ + 1) / 2));
    }

    RC_ExpectationCovarianceStatistic_2d(x, ix, y, iy, weights, subset, block, 
                                         varonly, REAL(L), REAL(E), REAL(V));
    UNPROTECT(1);
    return(ans);
}
