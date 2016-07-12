
void C_contrasts_marginal_maxabsstat(double *linstat, double *expect, double *covar,
                                     double *contrasts, int P, int Q, int Ncontrasts,
                                     double tol, int *wmax, double *maxstat);
void C_contrasts_marginal_quadform(double *linstat, double *expect, double *covar,
                                   double *contrasts, int P, int Q, int Ncontrasts,
                                   double tol, int *wmax, double *maxstat);
void C_ordered_maxabsstat_X
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
);

void C_ordered_quadform_X
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
);
