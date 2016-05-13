
void C_MPinv_sym (SEXP x, SEXP tol, double *dMP, int *rank);
int NLEVELS(SEXP x);
int NROW(SEXP x);
int NCOL(SEXP x);
void C_kronecker(const double *A, const int m, const int n,
                 const double *B, const int r, const int s, int overwrite,
                 double *ans);
void C_kronecker_sym(const double *A, const int m, 
                     const double *B, const int r, int overwrite,
                     double *ans);
