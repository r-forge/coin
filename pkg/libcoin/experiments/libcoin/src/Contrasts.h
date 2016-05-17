
void C_contrasts_marginal_maxabsstat(double *linstat, double *expect, double *covar,
                                     double *contrasts, int P, int Q, int Ncontrasts,
                                     double tol, int *wmax, double *maxstat);
void C_contrasts_marginal_quadform(double *linstat, double *expect, double *covar,
                                   double *contrasts, int P, int Q, int Ncontrasts,
                                   double tol, int *wmax, double *maxstat);
