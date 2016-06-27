
double C_quadform(int PQ, double *linstat, double *expect, double *MPinv);
double CR_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv);
double C_maxstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol);
double C_maxstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol);
double C_minstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol);
double C_minstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol);
double C_maxabsstat_Covariance(int PQ, double *linstat, double *expect, double *covar, double tol);
double C_maxabsstat_Variance(int PQ, double *linstat, double *expect, double *var, double tol);
double CR_maxstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol);
double CR_minstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol);
double CR_maxabsstat(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP tol);
