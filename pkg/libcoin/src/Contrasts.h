
void C_contrasts_marginal(double *linstat, double *expect, double *covar,
                          double *contrasts, int P, int Q, int teststat, int Ncontrasts,
                          double tol, int *wmax, double *maxstat);
void C_ordered_Xfactor(double *linstat, double *expect, double *covar, int P,
                       int Q, double *ExpX, int minbucket, double tol, int teststat,
                       int *wmax, double *maxstat);
void C_ordered_Xfactor_varonly(double *linstat, double *expect, double *varinf, int P,
                               int Q, double *ExpX, int minbucket, double tol, int teststat,
                               int *wmax, double *maxstat);
                               
                                                              