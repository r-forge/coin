
extern int C_get_P(SEXP LECV); 
extern int C_get_Q(SEXP LECV);
extern int C_get_Lb(SEXP LECV);
extern int C_get_varonly(SEXP LECV);
extern int C_get_Xfactor(SEXP LECV);
extern double* C_get_LinearStatistic(SEXP LECV);
extern double* C_get_Expectation(SEXP LECV);
extern double* C_get_Covariance(SEXP LECV);
extern double* C_get_MPinv(SEXP LECV);
extern double* C_get_Variance(SEXP LECV);
extern double* C_get_ExpectationX(SEXP LECV);
extern double* C_get_ExpectationInfluence(SEXP LECV);
extern double* C_get_CovarianceInfluence(SEXP LECV);
extern double* C_get_VarianceInfluence(SEXP LECV);
extern double* C_get_Work(SEXP LECV);
extern int* C_get_TableBlock(SEXP LECV);
extern int* C_get_Sumweights(SEXP LECV);
extern int* C_get_Table(SEXP LECV);
extern int* C_get_dimTable (SEXP LECV);
extern SEXP R_init_LECV(SEXP P, SEXP Q, SEXP varonly, SEXP Lb, SEXP Xfactor);
extern SEXP R_init_LECV_2d(SEXP P, SEXP Q, SEXP varonly, SEXP Lx, SEXP Ly, SEXP Lb, SEXP Xfactor);
