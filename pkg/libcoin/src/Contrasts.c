
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"

void C_contrasts_marginal_maxabsstat(double *linstat, double *expect, double *covar,
                                     double *contrasts, int P, int Q, int Ncontrasts,
                                     double tol, int *wmax, double *maxstat) {
                         
    double *mlinstat, *mexpect, *mvar, *mtmp, tmp;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mvar = Calloc(Q, double);
    mtmp = Calloc(P, double);
    wmax[0] = 0;
    maxstat[0] = 0.0;
    
    for (int i = 0; i < Ncontrasts; i++) {

        for (int q = 0; q < Q; q++) {
            mlinstat[q] = 0.0;
            mexpect[q] = 0.0;
            mvar[q] = 0.0;

            for (int p = 0; p < P; p++) {
                mlinstat[q] += contrasts[p + i * P] * linstat[q * P + p];
                mexpect[q] += contrasts[p + i * P] * expect[q * P + p];
            }
            /* only variances needed */
            for (int p = 0; p < P; p++) {
                mtmp[p] = 0.0;
                for (int pp = 0; pp < P; pp++)
                    mtmp[p] += contrasts[pp + i * P] * covar[S(pp + q * P, p + P * q, P * Q)];
            }
            for (int p = 0; p < P; p++)
                mvar[q] += contrasts[p + i * P] * mtmp[p];
        }

        tmp = C_maxabsstat_Variance(Q, mlinstat, mexpect, mvar, tol);
        if (tmp > maxstat[0]) {
            wmax[0] = i;
            maxstat[0] = tmp;
        }
    }
    Free(mlinstat); Free(mexpect); Free(mvar); Free(mtmp);
}

void C_contrasts_marginal_quadform(double *linstat, double *expect, double *covar,
                                   double *contrasts, int P, int Q, int Ncontrasts,
                                   double tol, int *wmax, double *maxstat) {
                         
    double *mlinstat, *mexpect, *mcovar, *mtmp, tmp, *MPinv;
    int rank, qq;
    
    mlinstat = Calloc(Q, double);
    mexpect = Calloc(Q, double);
    mcovar = Calloc(Q * (Q + 1) / 2, double);
    MPinv = Calloc(Q * (Q + 1) / 2, double);
    mtmp = Calloc(P, double);
    wmax[0] = 0;
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
                        mtmp[p] += contrasts[pp + i * P] * covar[S(pp + q * P, p + P * qq, P * Q)];
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
