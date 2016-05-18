
void C_MPinv_sym (double *x, int n, double tol, double *dMP, int *rank);
int NLEVELS(SEXP x);
int NROW(SEXP x);
int NCOL(SEXP x);
void C_kronecker(const double *A, const int m, const int n,
                 const double *B, const int r, const int s, int overwrite,
                 double *ans);
void C_kronecker_sym(const double *A, const int m, 
                     const double *B, const int r, int overwrite,
                     double *ans);
void
rcont2(int *nrow, int *ncol,
       int *nrowt, int *ncolt, int *ntotal,
              double *fact, int *jwork, int *matrix);
              