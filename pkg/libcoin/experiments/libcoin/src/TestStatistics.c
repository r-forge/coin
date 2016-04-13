
#include "libcoin.h"

double C_quadform(int PQ, double *linstat, double *expect, double *MPinv)
{
    int qPQ;
    double ans = 0.0, tmp = 0.0;
    
    for (int q = 0; q < PQ; q++) {
        qPQ = q * PQ;
        tmp = 0.0;
        for (int p = 0; p < PQ; p++)
            tmp += (linstat[p] - expect[p]) * MPinv[qPQ + p];
        ans += tmp * (linstat[q] - expect[q]);
    }
    return(ans);
}

double CR_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv)
{
    return(C_quadform(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                      REAL(Expectation), REAL(MPinv)));
}
