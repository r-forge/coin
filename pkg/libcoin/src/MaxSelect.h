
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
);

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
);

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
);

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
);

