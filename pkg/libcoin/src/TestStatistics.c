
#include "libcoin_internal.h"

double C_maxstand_Covariance(int PQ, double *linstat, double *expect, double *covar_sym, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_maxstand_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_minstand_Covariance(int PQ, double *linstat, double *expect, double *covar_sym, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

double C_minstand_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = (linstat[p] - expect[p]) / sqrt(var[p]);
        if (tmp < ans) ans = tmp;
    }
    return(ans);
}

double C_maxabsstand_Covariance(int PQ, double *linstat, double *expect, double *covar_sym, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (covar_sym[S(p, p, PQ)] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(covar_sym[S(p, p, PQ)]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_maxabsstand_Variance(int PQ, double *linstat, double *expect, double *var, double tol)
{

    double ans = 0.0, tmp = 0.0;
    
    for (int p = 0; p < PQ; p++) {
        tmp = 0.0;
        if (var[p] > tol)
            tmp = fabs((linstat[p] - expect[p]) / sqrt(var[p]));
        if (tmp > ans) ans = tmp;
    }
    return(ans);
}

double C_quadform(int PQ, double *linstat, double *expect, double *MPinv_sym)
{
    double ans = 0.0, tmp = 0.0;
    
    for (int q = 0; q < PQ; q++) {
        tmp = 0.0;
        for (int p = 0; p < PQ; p++)
            tmp += (linstat[p] - expect[p]) * MPinv_sym[S(p, q, PQ)];
        ans += tmp * (linstat[q] - expect[q]);
    }
    return(ans);
}

double C_maxtype(int PQ, double *linstat, double *expect, double *covar, int varonly, 
                 double tol, int alternative) 
{
    double ret;
    
    if (varonly) {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Variance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Variance(PQ, linstat, expect, covar, tol);
        } 
    } else {
        if (alternative ==  ALTERNATIVE_twosided) {
            ret = C_maxabsstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_less) {
            ret = C_minstand_Covariance(PQ, linstat, expect, covar, tol);
        } else if (alternative == ALTERNATIVE_greater) {
            ret = C_maxstand_Covariance(PQ, linstat, expect, covar, tol);
        } 
    }
    return(ret);
}

void C_standardise(int PQ, double *linstat, double *expect, double *covar, int varonly,
                   double tol)
{
    double var;
    
    for (int p = 0; p < PQ; p++) {
        if (varonly) {
            var = covar[p];
        } else {
            var = covar[S(p, p, PQ)];
        }
        if (var > tol) {
            linstat[p] = (linstat[p] - expect[p]) / sqrt(var);
        } else {
            linstat[p] = NAN;
        }
    }
}
