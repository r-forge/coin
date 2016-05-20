
#include "libcoin_internal.h"
#include "LinearStatistic.h"
#include "Utils.h"
#include "Sums.h"
#include "Tables.h"
#include "Distributions.h"
#include "MemoryAccess.h"

void RC_ExpectationCovarianceStatistic(SEXP x, SEXP y, SEXP weights, 
                                       SEXP subset, SEXP block, 
                                       SEXP ans)
{
    int N, P, Q, Nlevel, *sumweights, *table, *subset_tmp, tmp, chk;
    double *ExpInf, *CovInf, *work, *L, *E, *V;

    P = C_get_P(ans);
    Q = C_get_Q(ans);
    L = C_get_LinearStatistic(ans);
    E = C_get_Expectation(ans);

    N = NROW(x);
    Nlevel = 1;
    if (LENGTH(block) > 0)
        Nlevel = NLEVELS(block);

    ExpInf = Calloc(Nlevel * Q, double);
    if (C_get_varonly(ans)) {
        V = C_get_Variance(ans);
       CovInf = Calloc(Nlevel * Q, double);
       work = Calloc(3 * P + 1, double);
    } else {
        V = C_get_Covariance(ans);
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

    if (C_get_varonly(ans)) {
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
    SEXP ans, P, Q; 
    
    PROTECT(P = ScalarInteger(0));
    PROTECT(Q = ScalarInteger(0));

    if (isInteger(x)) {
        INTEGER(P)[0] = NLEVELS(x);
    } else {
        INTEGER(P)[0] = NCOL(x);
    }
    INTEGER(Q)[0] = NCOL(y);
    
    PROTECT(ans = R_init_LECV(P, Q, varonly));

    RC_ExpectationCovarianceStatistic(x, y, weights, subset, block, ans);

    UNPROTECT(3);
    return(ans);
}



void RC_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy,
                                          SEXP weights, SEXP subset, SEXP block, 
                                          SEXP ans)
{

    int N, P, Q, Lxp1, Lyp1, Lb, *btab, *csum, *rsum, *table, *table2d, sw, *iix;
    double *L, *E, *V, *ExpInf, *CovInf, *ExpX, *CovX, *work;

    N = LENGTH(ix);
    P = C_get_P(ans);
    Q = C_get_Q(ans);

    L = C_get_LinearStatistic(ans);
    E = C_get_Expectation(ans);
    table = C_get_Table(ans);

    Lxp1 = C_get_dimTable(ans)[0];
    Lyp1 = C_get_dimTable(ans)[1];
    Lb = C_get_dimTable(ans)[2];

    table2d = Calloc(Lxp1 * Lyp1, int);
    csum = Calloc(Lyp1, int);
    rsum = Calloc(Lxp1, int);
    
    for (int i = 0; i < Lxp1 * Lyp1; i++)
        table2d[i] = 0;
    for (int b = 0; b < Lb; b++) {
        for (int i = 0; i < Lxp1; i++) {
            for (int j = 0; j < Lyp1; j++)
                table2d[j * Lxp1 + i] += table[b * Lxp1 * Lyp1 + j * Lxp1 + i];
        }
    }


    ExpInf = Calloc(Q, double);
    ExpX = Calloc(P, double);
    if (C_get_varonly(ans)) {
        V = C_get_Variance(ans);
        CovInf = Calloc(Q, double);
        CovX = Calloc(P, double);
        work = Calloc(P, double);
    } else {
        V = C_get_Covariance(ans);
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

    for (int b = 0; b < Lb; b++) {
        btab = table + Lxp1 * Lyp1 * b;
        C_colSums_i(btab, Lxp1, Lyp1, csum); 
        csum[0] = 0; /* NA */   
        C_rowSums_i(btab, Lxp1, Lyp1, rsum);    
        rsum[0] = 0; /* NA */
        sw = 0;
        for (int i = 1; i < Lxp1; i++) sw += rsum[i];
        C_ExpectationInfluence_weights(REAL(y), NROW(y), Q, csum, sw, ExpInf);
        if (LENGTH(x) == 0) {
            for (int p = 0; p < P; p++)
                ExpX[p] = (double) rsum[p + 1];
        } else {
            C_ExpectationX_weights(REAL(x), NROW(x), P, rsum, ExpX);
        }
        C_ExpectationLinearStatistic(P, Q, ExpInf, ExpX, b, E);
        if (C_get_varonly(ans)) {
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
    Free(table2d); Free(csum); Free(rsum); Free(ExpInf); Free(ExpX);
    Free(CovInf); Free(CovX); Free(work);
}        


SEXP R_ExpectationCovarianceStatistic_2d(SEXP x, SEXP ix, SEXP y, SEXP iy, 
                                         SEXP weights, SEXP subset, SEXP block, 
                                         SEXP varonly)
{
    SEXP ans, P, Q, Lx, Ly, Lb;

    PROTECT(P = ScalarInteger(0));
    PROTECT(Q = ScalarInteger(0));
    PROTECT(Lx = ScalarInteger(0));
    PROTECT(Ly = ScalarInteger(0));
    PROTECT(Lb = ScalarInteger(0));
    
    if (LENGTH(x) == 0) {
        INTEGER(P)[0] = NLEVELS(ix);
    } else {
        INTEGER(P)[0] = NCOL(x);
    }
    INTEGER(Q)[0] = NCOL(y);

    INTEGER(Lb)[0] = 1;
    if (LENGTH(block) > 0)
        INTEGER(Lb)[0] = NLEVELS(block);
        
    INTEGER(Lx)[0] = NLEVELS(ix);
    INTEGER(Ly)[0] = NLEVELS(iy);

    PROTECT(ans = R_init_LECV_2d(P, Q, varonly, Lx, Ly, Lb));

    RC_2dtable(ix, iy, weights, subset, block, C_get_Table(ans));
    RC_ExpectationCovarianceStatistic_2d(x, ix, y, iy, weights, 
                                         subset, block, ans);

    UNPROTECT(6);
    return(ans);
}

SEXP R_PermutedLinearStatistic_2d(SEXP LEV, SEXP x, SEXP ix, SEXP y, SEXP iy, 
                                  SEXP block, SEXP B) {

    SEXP ans;
    int N, P, Q, PQ, Lb, Lx, Ly, *csum, *rsum, *ntotal, *table, *jwork, *rtable, *rtable2, maxn = 0, Lxp1, Lyp1;
    double *fact, *linstat, *blinstat;
    
    N = LENGTH(ix);
    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
    Lxp1 = C_get_dimTable(LEV)[0];
    Lyp1 = C_get_dimTable(LEV)[1];
    Lx = Lxp1 - 1;
    Ly = Lyp1 - 1;
    Lb = C_get_dimTable(LEV)[2];

    table = C_get_Table(LEV);

    PROTECT(ans = allocMatrix(REALSXP, PQ, INTEGER(B)[0]));
    
    csum = Calloc(Lyp1 * Lb, int);
    rsum = Calloc(Lxp1 * Lb, int);
    ntotal = Calloc(Lb, int);
    rtable = Calloc(Lxp1 * Lyp1, int);
    rtable2 = Calloc(NLEVELS(ix) * NLEVELS(iy) , int);
    linstat = Calloc(PQ, double);
    jwork = Calloc(Lyp1, int);

    for (int b = 0; b < Lb; b++) {
        C_colSums_i(table + Lxp1 * Lyp1 * b, 
                    NLEVELS(ix) + 1, NLEVELS(iy) + 1, csum + Lyp1 * b); 
        csum[Lyp1 * b] = 0; /* NA */   
        C_rowSums_i(table + Lxp1 * Lyp1 * b, 
                    Lxp1, Lyp1, rsum + Lxp1 * b);
        rsum[Lxp1 * b] = 0; /* NA */
        ntotal[b] = 0;
        for (int i = 1; i < Lxp1; i++) 
            ntotal[b] += rsum[Lxp1 * b + i];
        if (ntotal[b] > maxn) maxn = ntotal[b];
    }
    
    fact = Calloc(maxn + 1, double);    
    /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
    fact[0] = fact[1] = 0.;
    for(int j = 2; j <= maxn; j++)
        fact[j] = fact[j - 1] + log(j);

    GetRNGstate(); 

    for (int i = 0; i < INTEGER(B)[0]; i++) {
        blinstat = REAL(ans) + PQ * i;
        for (int p = 0; p < PQ; p++) {
            blinstat[p] = 0;
            linstat[p] = 0;
        }
        for (int p = 0; p < Lxp1 * Lyp1; p++)
            rtable[p] = 0;
            
        for (int b = 0; b < Lb; b++) {
            
            rcont2(&Lx, &Ly, rsum + Lxp1 * b + 1, 
                   csum + Lyp1 * b + 1, ntotal + b, fact, jwork, rtable2);

            for (int j1 = 1; j1 <= NLEVELS(ix); j1++) {
                for (int j2 = 1; j2 <= NLEVELS(iy); j2++) 
                    rtable[j2 * Lxp1 + j1] = rtable2[(j2 - 1) * NLEVELS(ix) + (j1 - 1)];
            }
            C_LinearStatistic_2d(x, Lxp1, P, REAL(y), NROW(y), Q, rtable, linstat);
            for (int p = 0; p < PQ; p++)
                blinstat[p] += linstat[p];
        }
    }
    
    PutRNGstate();
    
    Free(csum); Free(rsum); Free(ntotal); Free(rtable); Free(rtable2); Free(linstat); Free(jwork); Free(fact);
    UNPROTECT(1);
    return(ans);
}
                  