
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "LinearStatistic.h"

void C_contrasts_marginal_maxabsstand
(
    double *linstat, 
    double *expect, 
    double *covar,
    double *contrasts, 
    int P, 
    int Q, 
    int Ncontrasts,
    double tol, 
    int *wmax, 
    double *maxstat
) {
                         
    double *mlinstat, *mexpect, *mvar, *mtmp, tmp;
    int PQ = P * Q, qPp;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mvar = Calloc(Q, double);
    mtmp = Calloc(P, double);
    wmax[0] = -1;
    maxstat[0] = 0.0;
    
    for (int i = 0; i < Ncontrasts; i++) {

        for (int q = 0; q < Q; q++) {
            mlinstat[q] = 0.0;
            mexpect[q] = 0.0;
            mvar[q] = 0.0;

            for (int p = 0; p < P; p++) {
                qPp = q * P + p;
                mlinstat[q] += contrasts[p + i * P] * linstat[qPp];
                mexpect[q] += contrasts[p + i * P] * expect[qPp];
                /* only variances needed */
                mtmp[p] = 0.0;
                for (int pp = 0; pp < P; pp++)
                    mtmp[p] += contrasts[pp + i * P] * 
                               covar[S(pp + q * P, qPp, PQ)];
            }
            for (int p = 0; p < P; p++)
                mvar[q] += contrasts[p + i * P] * mtmp[p];
        }

        tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                        ALTERNATIVE_twosided);
        
        if (tmp > maxstat[0]) {
            wmax[0] = i;
            maxstat[0] = tmp;
        }
    }
    Free(mlinstat); Free(mexpect); Free(mvar); Free(mtmp);
}

void C_contrasts_marginal_quadform
(
    double *linstat, 
    double *expect, 
    double *covar,
    double *contrasts, 
    int P, 
    int Q, 
    int Ncontrasts,
    double tol, 
    int *wmax, 
    double *maxstat
) {
                         
    double *mlinstat, *mexpect, *mcovar, *mtmp, tmp, *MPinv;
    int rank, qq;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mcovar = Calloc(Q * (Q + 1) / 2, double);
    MPinv = Calloc(Q * (Q + 1) / 2, double);
    mtmp = Calloc(P, double);
    wmax[0] = -1;
    maxstat[0] = 0.0;
    
    for (int i = 0; i < Ncontrasts; i++) {

        for (int q = 0; q < Q; q++) {
            mlinstat[q] = 0.0;
            mexpect[q] = 0.0;
            for (int qq = 0; qq <= q; qq++)
                mcovar[S(q, qq, Q)] = 0.0;

            for (int p = 0; p < P; p++) {
                mlinstat[q] += contrasts[p + i * P] * linstat[q * P + p];
                mexpect[q] += contrasts[p + i * P] * expect[q * P + p];
            }

            for (qq = 0; qq <= q; qq++) {
                for (int p = 0; p < P; p++) {
                    mtmp[p] = 0.0;
                    for (int pp = 0; pp < P; pp++)
                        mtmp[p] += contrasts[pp + i * P] * 
                                   covar[S(pp + q * P, p + P * qq, P * Q)];
                }
                for (int p = 0; p < P; p++)
                    mcovar[S(q, qq, Q)] += contrasts[p + i * P] * mtmp[p];
            }
        }

        C_MPinv_sym(mcovar, Q, tol, MPinv, &rank);
        tmp = C_quadform(Q, mlinstat, mexpect, MPinv);
        
        if (tmp > maxstat[0]) {
            wmax[0] = i;
            maxstat[0] = tmp;
        }
    }
    Free(mlinstat); Free(mexpect); Free(mcovar); Free(mtmp);
    Free(MPinv);
}

void C_contrasts_marginal
(
    double *linstat, 
    double *expect, 
    double *covar,
    double *contrasts, 
    int P, 
    int Q, 
    int teststat, 
    int Ncontrasts,
    double tol, 
    int *wmax, 
    double *maxstat
) {

    if (teststat == TESTSTAT_maxtype) {
        C_contrasts_marginal_maxabsstand(linstat, expect, covar,
                                         contrasts, P, Q, Ncontrasts,
                                         tol, wmax, maxstat);
    } else {
        C_contrasts_marginal_quadform(linstat, expect, covar,
                                      contrasts, P, Q, Ncontrasts,
                                      tol, wmax, maxstat);
    }
}

void C_ordered_maxabsstand_Xfactor
(
    double *linstat, 
    double *expect, 
    double *covar, 
    int P, 
    int Q,
    double *ExpX, 
    int minbucket, 
    double tol, 
    int *wmax, 
    double *maxstat
) {

    double *mlinstat, *mexpect, *mvar, tmp, sumleft, sumright;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mvar = Calloc(Q, double);
    wmax[0] = -1;
    maxstat[0] = 0.0;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        mvar[q] = 0.0;
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
            mexpect[q] += expect[q * P + p];
            /* only variances needed */
            for (int pp = 0; pp < p; pp++)
                mvar[q] += 2 * covar[S(pp + q * P, p + P * q, P * Q)];
            mvar[q] += covar[S(p + q * P, p + P * q, P * Q)];
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket) && 
            (ExpX[p] > 0)) {

            tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                            ALTERNATIVE_twosided);
            
            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }
        }
    }
    Free(mlinstat); Free(mexpect); Free(mvar);    
}

void C_ordered_quadform_Xfactor
(
    double *linstat, 
    double *expect, 
    double *covar, 
    int P, 
    int Q,
    double *ExpX, 
    int minbucket, 
    double tol, 
    int *wmax, 
    double *maxstat
) {

    double *mlinstat, *mexpect, *mcovar, *MPinv, tmp, sumleft, sumright;
    int rank;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mcovar = Calloc(Q * (Q + 1) / 2, double);
    MPinv = Calloc(Q * (Q + 1) / 2, double);
    wmax[0] = -1;
    maxstat[0] = 0.0;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        for (int qq = 0; qq <= q; qq++)
            mcovar[S(q, qq, Q)] = 0.0;
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
            mexpect[q] += expect[q * P + p];
            for (int qq = 0; qq <= q; qq++) {
                for (int pp = 0; pp < p; pp++)
                    mcovar[S(q, qq, Q)] += 
                        2 * covar[S(pp + q * P, p + P * qq, P * Q)];
                mcovar[S(q, qq, Q)] += 
                    covar[S(p + q * P, p + P * qq, P * Q)];
            }
        }
                 
        if ((sumleft >= minbucket) && 
            (sumright >= minbucket) && 
            (ExpX[p] > 0)) {
        
            C_MPinv_sym(mcovar, Q, tol, MPinv, &rank);
            tmp = C_quadform(Q, mlinstat, mexpect, MPinv);

            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }
        }
    }
    Free(mlinstat); Free(mexpect); Free(mcovar); Free(MPinv);   
}


void C_ordered_Xfactor
(
    double *linstat, 
    double *expect, 
    double *covar, 
    int P, 
    int Q, 
    double *ExpX, 
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat
) {

    if (teststat == TESTSTAT_maxtype) {
        C_ordered_maxabsstand_Xfactor(linstat, expect, covar, P, Q, ExpX, 
                                      minbucket, tol, wmax, maxstat);
    } else {
        C_ordered_quadform_Xfactor(linstat, expect, covar, P, Q, ExpX, 
                                   minbucket, tol, wmax, maxstat);
    }
}

void C_ordered_Xfactor_varonly
(
    double *linstat, 
    double *expect, 
    double *varinf, 
    int P,
    int Q, 
    double *ExpX, 
    int minbucket, 
    double tol, 
    int teststat,
    int *wmax, 
    double *maxstat
) {

    double *mlinstat, *mexpect, *mvar, tmp, sumleft, sumright, Ptmp;
    int sw;

    /* quadform needs covinf not varinf for Q > 1*/
    if (teststat != TESTSTAT_maxtype)
        error("only maxtype implemented");
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mvar = Calloc(Q, double);
    wmax[0] = -1;
    maxstat[0] = 0.0;

    for (int q = 0; q < Q; q++) {
        mlinstat[q] = 0.0;
        mexpect[q] = 0.0;
        mvar[q] = 0.0;
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
            mexpect[q] += expect[q * P + p];
            /* does not work with blocks! */
            C_VarianceLinearStatistic(1, Q, varinf, &sumleft, &sumleft,
                                      sw, &Ptmp, 0, mvar);
        }

        if ((sumleft >= minbucket) && 
            (sumright >= minbucket) && 
            (ExpX[p] > 0)) {

            tmp = C_maxtype(Q, mlinstat, mexpect, mvar, 1, tol, 
                            ALTERNATIVE_twosided);
            
            if (tmp > maxstat[0]) {
                wmax[0] = p;
                maxstat[0] = tmp;
            }
        }
    }
    Free(mlinstat); Free(mexpect); Free(mvar);    
}
