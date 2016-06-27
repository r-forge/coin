
int C_get_P(SEXP LECV); 
int C_get_Q(SEXP LECV);
int C_get_varonly(SEXP LECV);
double* C_get_LinearStatistic(SEXP LECV);
double* C_get_Expectation(SEXP LECV);
double* C_get_Covariance(SEXP LECV);
double* C_get_MPinv(SEXP LECV);
double* C_get_Variance(SEXP LECV);
double* C_get_ExpectationX(SEXP LECV);
double* C_get_ExpectationInfluence(SEXP LECV);
double* C_get_CovarianceInfluence(SEXP LECV);
double* C_get_Work(SEXP LECV);
int* C_get_TableBlock(SEXP LECV);
int* C_get_Sumweights(SEXP LECV);
int* C_get_Table(SEXP LECV);
int* C_get_dimTable (SEXP LECV);
SEXP R_init_LECV(SEXP P, SEXP Q, SEXP varonly, SEXP Lb);
SEXP R_init_LECV_2d(SEXP P, SEXP Q, SEXP varonly, SEXP Lx, SEXP Ly, SEXP Lb);