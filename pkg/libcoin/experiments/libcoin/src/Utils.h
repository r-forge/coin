
void CR_La_svd(int dim, char *jobu, char *jobv, double *x, 
               double *s, double *u, double *v);
void CR_svd (SEXP x, SEXP s, SEXP u, SEXP v);
void C_MPinv (SEXP x, double tol, SEXP s, SEXP u, SEXP v, double *MPinv, 
              int *rank);
