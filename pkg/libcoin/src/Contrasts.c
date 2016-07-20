
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "LinearStatistic.h"
#include "Distributions.h"

void C_ordered_Xfactor_block
(
    double *linstat, 
    double *expect, 
    double *covar,
    int P,
    int Q, 
    double *ExpX, 
    int B,
    double* blinstat,
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat,
    double *pval,
    int lower,
    int give_log
) {

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, *bmaxstat, 
           tmp, sumleft, sumright;
    int rank, PQ = P * Q, greater;

    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    if (teststat == TESTSTAT_maxtype) {
       mvar = Calloc(Q, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
    }
    if (B > 0) {
        mblinstat = Calloc(Q * B, double);
        bmaxstat = Calloc(B, double);
    }
       
    wmax[0] = NA_INTEGER;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        if (teststat == TESTSTAT_maxtype)
            mvar[q] = 0.0;
        for (int b = 0; b < B; b++) mblinstat[q + b * Q] = 0.0;
    }
    if (teststat == TESTSTAT_quadform) {
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            mcovar[q] = 0.0;
    }
    
    sumleft = 0.0;                        
    sumright = 0.0;
    for (int p = 0; p < P; p++) 
        sumright += ExpX[p];
                 
    for (int p = 0; p < P; p++) {
        sumleft += ExpX[p];
        sumright -= ExpX[p];

        for (int q = 0; q < Q; q++) {
            mlinstat[q] += linstat[q * P + p];
            for (int b = 0; b < B; b++)
                mblinstat[q + b * Q] += blinstat[q * P + p + b * PQ];
            mexpect[q] += expect[q * P + p];
            if (teststat == TESTSTAT_maxtype) {
                for (int pp = 0; pp < p; pp++)
                    mvar[q] += 2 * covar[S(pp + q * P, p + P * q, P * Q)];
                mvar[q] += covar[S(p + q * P, p + P * q, P * Q)];
            } else {
                for (int qq = 0; qq <= q; qq++) {
                    for (int pp = 0; pp < p; pp++)
                        mcovar[S(q, qq, Q)] += 
                            2 * covar[S(pp + q * P, p + P * qq, P * Q)];
                    mcovar[S(q, qq, Q)] += 
                        covar[S(p + q * P, p + P * qq, P * Q)];
                }
            }
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket) && 
            (ExpX[p] > 0)) {

            if (teststat == TESTSTAT_maxtype) {
                tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, mlinstat, mexpect, mcovar);
            }

            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }

            for (int b = 0; b < B; b++) {
                if (teststat == TESTSTAT_maxtype) {
                    tmp = C_maxtype(Q, mblinstat + b * Q, mexpect, mvar, 1, tol, 
                                    ALTERNATIVE_twosided);
                } else {
                    tmp = C_quadform(Q, mblinstat + b * Q, mexpect, mcovar);
                }
                if (tmp > bmaxstat[b])
                    bmaxstat[b] = tmp;
            }
        }
        if (B > 0) {
            greater = 0;
            for (int b = 0; b < B; b++) {
                if (bmaxstat[b] > maxstat[0]) greater++;
            }
            pval[0] = C_perm_pvalue(greater, B, lower, give_log);
        }
    }
    Free(mlinstat); Free(mexpect); 
    if (B > 0) {
        Free(mblinstat); Free(bmaxstat);
    }
    if (teststat == TESTSTAT_maxtype) {
        Free(mvar);
    } else {
        Free(mcovar); Free(mMPinv);
    }
}

void C_ordered_Xfactor
(
    double *linstat, 
    double *expect, 
    double *varinf,
    double *covinf,
    int P,
    int Q, 
    double *ExpX, 
    int B,
    double *blinstat,
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat,
    double *pval,
    int lower,
    int give_log
) {

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, *bmaxstat, 
           tmp, sumleft, sumright, Ptmp;
    int sw, rank, PQ = P * Q, greater;

    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    if (teststat == TESTSTAT_maxtype) {
       mvar = Calloc(Q, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
    }
    if (B > 0) {
        mblinstat = Calloc(Q * B, double);
        bmaxstat = Calloc(B, double);
    }
       
    wmax[0] = NA_INTEGER;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        if (teststat == TESTSTAT_maxtype)
            mvar[q] = 0.0;
        for (int b = 0; b < B; b++) mblinstat[q + b * Q] = 0.0;
    }
    if (teststat == TESTSTAT_quadform) {
        for (int q = 0; q < Q * (Q + 1) / 2; q++)
            mcovar[q] = 0.0;
    }
    
    sumleft = 0.0;                        
    sumright = 0.0;
    for (int p = 0; p < P; p++) 
        sumright += ExpX[p];
    sw = sumright;
                 
    for (int p = 0; p < P; p++) {
        sumleft += ExpX[p];
        sumright -= ExpX[p];

        for (int q = 0; q < Q; q++) {
            mlinstat[q] += linstat[q * P + p];
            for (int b = 0; b < B; b++)
                mblinstat[q + b * Q] += blinstat[q * P + p + b * PQ];
            mexpect[q] += expect[q * P + p];
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket) && 
            (ExpX[p] > 0)) {

            /* does not work with blocks! */
            if (teststat == TESTSTAT_maxtype) {
                C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                                          sw, &Ptmp, 0, mvar);
            } else {
                C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                            sw, &Ptmp, 0, mcovar);
            }

            if (teststat == TESTSTAT_maxtype) {
                tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, mlinstat, mexpect, mcovar);
            }

            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }

            for (int b = 0; b < B; b++) {
                if (teststat == TESTSTAT_maxtype) {
                    tmp = C_maxtype(Q, mblinstat + b * Q, mexpect, mvar, 1, tol, 
                                    ALTERNATIVE_twosided);
                } else {
                    tmp = C_quadform(Q, mblinstat + b * Q, mexpect, mcovar);
                }
                if (tmp > bmaxstat[b])
                    bmaxstat[b] = tmp;
            }
        }
        if (B > 0) {
            greater = 0;
            for (int b = 0; b < B; b++) {
                if (bmaxstat[b] > maxstat[0]) greater++;
            }
            pval[0] = C_perm_pvalue(greater, B, lower, give_log);
        }
    }
    Free(mlinstat); Free(mexpect); 
    if (B > 0) {
        Free(mblinstat); Free(bmaxstat);
    }
    if (teststat == TESTSTAT_maxtype) {
        Free(mvar);
    } else {
        Free(mcovar); Free(mMPinv);
    }
}

void C_unordered_Xfactor_block
(
    double *linstat, 
    double *expect, 
    double *covar,
    int P,
    int Q, 
    double *ExpX, 
    int B,
    double* blinstat,
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat,
    double *pval,
    int lower,
    int give_log
) {

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, *bmaxstat, 
           tmp, sumleft, sumright, *mtmp;
    int rank, PQ = P * Q, qPp, greater, nc, *levels, Pnonzero, *indl, *contrast;

    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mtmp = Calloc(P, double);
    if (teststat == TESTSTAT_maxtype) {
       mvar = Calloc(Q, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
    }
    if (B > 0) {
        mblinstat = Calloc(Q * B, double);
        bmaxstat = Calloc(B, double);
    }

    for (int p = 0; p < P; p++) wmax[p] = NA_INTEGER;       

    contrast = Calloc(P, int);
    for (int p = 0; p < P; p++) {
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

    int mi = 1;  
    for (int l = 1; l < Pnonzero; l++) mi *= 2;
    indl = Calloc(Pnonzero, int);
    for (int p = 0; p < Pnonzero; p++) indl[p] = 0;
       
    for (int j = 1; j < mi; j++) { /* go though all splits */
    
        /* indl determines if level p is left or right */
        int jj = j;
        for (int l = 1; l < Pnonzero; l++) {
            indl[l] = (jj%2);
            jj /= 2;
        }

        sumleft = 0.0;
        sumright = 0.0;
        for (int p = 0; p < P; p++) contrast[p] = 0;
        for (int p = 0; p < Pnonzero; p++) {
            sumleft += indl[p] * ExpX[levels[p]];
            sumright += (1 - indl[p]) * ExpX[levels[p]];
            contrast[levels[p]] = indl[p];
        }
                                                                                                                        
        if (teststat == TESTSTAT_maxtype) {
            for (int q = 0; q < Q; q++) {
                mlinstat[q] = 0.0;
                mexpect[q] = 0.0;
                mvar[q] = 0.0;
                for (int b = 0; b < B; b++)   
                    mblinstat[q + b * Q] = 0.0;
                
                for (int p = 0; p < P; p++) {
                    qPp = q * P + p;
                    mlinstat[q] += contrast[p] * linstat[qPp];
                    mexpect[q] += contrast[p] * expect[qPp];
                    for (int b = 0; b < B; b++)   
                        mblinstat[q + b * Q] += contrast[p] * blinstat[q * P + p + b * PQ];
                    mtmp[p] = 0.0;
                    for (int pp = 0; pp < P; pp++)
                        mtmp[p] += contrast[pp] * 
                                   covar[S(pp + q * P, qPp, PQ)];
                }
                for (int p = 0; p < P; p++)
                    mvar[q] += contrast[p] * mtmp[p];
            }
        } else {
            for (int q = 0; q < Q; q++) {
                mlinstat[q] = 0.0;
                mexpect[q] = 0.0;
                for (int qq = 0; qq <= q; qq++)
                    mcovar[S(q, qq, Q)] = 0.0;

                for (int p = 0; p < P; p++) {
                    qPp = q * P + p;
                    mlinstat[q] += contrast[p] * linstat[qPp];
                    mexpect[q] += contrast[p] * expect[qPp];
                }
                for (int qq = 0; qq <= q; qq++) {
                    for (int p = 0; p < P; p++) {
                        mtmp[p] = 0.0;
                        for (int pp = 0; pp < P; pp++)
                            mtmp[p] += contrast[pp] * 
                                       covar[S(pp + q * P, p + P * qq, P * Q)];
                    }
                    for (int p = 0; p < P; p++)
                        mcovar[S(q, qq, Q)] += contrast[p] * mtmp[p];
                }
            }
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket)) {

            if (teststat == TESTSTAT_maxtype) {
                tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, mlinstat, mexpect, mcovar);
            }

            if (tmp > maxstat[0]) {
                for (int p = 0; p < Pnonzero; p++)
                    wmax[levels[p]] = contrast[p];
                maxstat[0] = tmp;
            }

            for (int b = 0; b < B; b++) {
                if (teststat == TESTSTAT_maxtype) {
                    tmp = C_maxtype(Q, mblinstat + b * Q, mexpect, mvar, 1, tol, 
                                    ALTERNATIVE_twosided);
                } else {
                    tmp = C_quadform(Q, mblinstat + b * Q, mexpect, mcovar);
                }
                if (tmp > bmaxstat[b])
                    bmaxstat[b] = tmp;
            }
        }
        if (B > 0) {
            greater = 0;
            for (int b = 0; b < B; b++) {
                if (bmaxstat[b] > maxstat[0]) greater++;
            }
            pval[0] = C_perm_pvalue(greater, B, lower, give_log);
        }
    }
    Free(mlinstat); Free(mexpect); Free(levels); Free(contrast); Free(indl); Free(mtmp);
    if (B > 0) {
        Free(mblinstat); Free(bmaxstat);
    }
    if (teststat == TESTSTAT_maxtype) {
        Free(mvar);
    } else {
        Free(mcovar); Free(mMPinv);
    }
}

void C_unordered_Xfactor
(
    double *linstat, 
    double *expect, 
    double *varinf,
    double *covinf,
    int P,
    int Q, 
    double *ExpX, 
    int B,
    double *blinstat,
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat,
    double *pval,
    int lower,
    int give_log
) {

    double *mlinstat, *mblinstat, *mexpect, *mvar, *mcovar, *mMPinv, *bmaxstat, 
           tmp, sumleft, sumright, Ptmp;
    int sw, rank, PQ = P * Q, qPp, greater, *levels, *indl, *contrast;

    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    if (teststat == TESTSTAT_maxtype) {
       mvar = Calloc(Q, double);
    } else {
       mcovar = Calloc(Q * (Q + 1) / 2, double);
       mMPinv = Calloc(Q * (Q + 1) / 2, double);
    }
    if (B > 0) {
        mblinstat = Calloc(Q * B, double);
        bmaxstat = Calloc(B, double);
    }

    for (int p = 0; p < P; p++) wmax[p] = NA_INTEGER;

    contrast = Calloc(P, int);
    int Pnonzero = 0;
    sw = 0;
    for (int p = 0; p < P; p++) {
        sw += ExpX[p];
        if (ExpX[p] > 0) Pnonzero++;
    }
    levels = Calloc(Pnonzero, int);
    int nc = 0;
    for (int p = 0; p < P; p++) {
        if (ExpX[p] > 0) {
            levels[nc] = p;
            nc++;
        }
    }
    
    if (Pnonzero >= 31)
        error("cannot search for unordered splits in >= 31 levels");

    for (int p = 0; p < P; p++) {
        if (ExpX[p] == 0.0) wmax[p] = NA_INTEGER;
    }

    int mi = 1;  
    for (int l = 1; l < Pnonzero; l++) mi *= 2;
    indl = Calloc(Pnonzero, int);
    for (int p = 0; p < Pnonzero; p++) indl[p] = 0;
       
    for (int j = 1; j < mi; j++) { /* go though all splits */
    
        /* indl determines if level p is left or right */
        int jj = j;
        for (int l = 1; l < Pnonzero; l++) {
            indl[l] = (jj%2);
            jj /= 2;
        }

        sumleft = 0.0;
        sumright = 0.0;
        for (int p = 0; p < Pnonzero; p++) {
            sumleft += indl[p] * ExpX[levels[p]];
            sumright += (1 - indl[p]) * ExpX[levels[p]];
            contrast[levels[p]] = indl[p];
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket)) {

            for (int q = 0; q < Q; q++) {
                mlinstat[q] = 0.0;
                mexpect[q] = 0.0;
                for (int b = 0; b < B; b++)   
                    mblinstat[q + b * Q] = 0.0;
                for (int p = 0; p < P; p++) {
                    qPp = q * P + p;
                    mlinstat[q] += contrast[p] * linstat[qPp];
                    mexpect[q] += contrast[p] * expect[qPp];
                    for (int b = 0; b < B; b++)   
                        mblinstat[q + b * Q] += contrast[p] * blinstat[q * P + p + b * PQ];

                }
            }

            /* does not work with blocks! */
            if (teststat == TESTSTAT_maxtype) {
                C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                                          sw, &Ptmp, 0, mvar);
            } else {
                C_CovarianceLinearStatistic(1, Q, covinf, &sumleft, &sumleft,
                                            sw, &Ptmp, 0, mcovar);
            }

            if (teststat == TESTSTAT_maxtype) {
                tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                                ALTERNATIVE_twosided);
            } else {
                C_MPinv_sym(mcovar, Q, tol, mMPinv, &rank);
                tmp = C_quadform(Q, mlinstat, mexpect, mcovar);
            }

            if (tmp > maxstat[0]) {
                for (int p = 0; p < Pnonzero; p++)
                    wmax[levels[p]] = contrast[p];
                maxstat[0] = tmp;
            }

            for (int b = 0; b < B; b++) {
                if (teststat == TESTSTAT_maxtype) {
                    tmp = C_maxtype(Q, mblinstat + b * Q, mexpect, mvar, 1, tol, 
                                    ALTERNATIVE_twosided);
                } else {
                    tmp = C_quadform(Q, mblinstat + b * Q, mexpect, mcovar);
                }
                if (tmp > bmaxstat[b])
                    bmaxstat[b] = tmp;
            }
        }
        if (B > 0) {
            greater = 0;
            for (int b = 0; b < B; b++) {
                if (bmaxstat[b] > maxstat[0]) greater++;
            }
            pval[0] = C_perm_pvalue(greater, B, lower, give_log);
        }
    }
    Free(mlinstat); Free(mexpect); Free(levels); Free(contrast); Free(indl);
    if (B > 0) {
        Free(mblinstat); Free(bmaxstat);
    }
    if (teststat == TESTSTAT_maxtype) {
        Free(mvar);
    } else {
        Free(mcovar); Free(mMPinv);
    }
}
