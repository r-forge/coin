
double C_quadform(int PQ, double *linstat, double *expect, double *MPinv);
double C_maxtype(int PQ, double *linstat, double *expect, double *covar, int varonly,
                 double tol, int alternative);
void C_standardise(int PQ, double *linstat, double *expect, double *covar, int varonly,
                   double tol);
                   