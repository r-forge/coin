
void C_LinearStatistic(SEXP x, int N, int P, double* y, int Q, 
                       int *weights, int *sumweights, 
                       int *subset, int *Nsubset, int Nlevel, 
                       double *PQ_ans);
void C_PermutedLinearStatistic(SEXP x, int N, int P, double* y, int Q,
                               int *perm, int *original, int Nperm,
                               double *PQ_ans);
void C_LinearStatistic_2d(SEXP x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans); 
void C_LinearStatisticXfactor_2d(int N, double *y, int M, int Q, 
                                 int *weights, double *Nm1Q_ans);
void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *ExpX,
                                  int add, double *PQ_ans);
void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *ExpX, 
                                 double *CovX, int sumweights, double *PP_tmp, int add,
                                 double *PQPQ_ans); 
void C_VarianceLinearStatistic(int P, int Q, double *VarInf, double *ExpX, 
                               double *VarX, int sumweights, double *P_tmp, int add,
                               double *PQ_ans); 
void C_ExpectationCoVarianceInfluence(double* y, int N, int Q,    
                                      int *weights, int *sumweights,
                                      int *subset, int *Nsubset, int Nlevel, int varonly,
                                      double *NlevelQ_ans, double *NlevelQQ_sym_ans);
void C_ExpectationCovarianceLinearStatistic(SEXP x, int N, int P, int Q,
                                            int *weights, int *sumweights,
                                            int *subset, int *Nsubset, int Nlevel, double *ExpXtotal,
                                            double *ExpInf, double *CovInf, double *work, double *PQ_ans,
                                            double *PQPQ_sym_ans);
void C_ExpectationVarianceLinearStatistic(SEXP x, int N, int P, int Q,
                                          int *weights, int *sumweights,
                                          int *subset, int *Nsubset, int Nlevel, double *ExpXtotal,
                                          double *ExpInf, double *VarInf, double *work, double *PQ_ans_Exp, 
                                          double *PQ_ans_Var);
                                          
void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights,
                                    int sumweights, double *Q_ans);
void C_ExpectationInfluence_weights_subset(double* y, int N, int Q,
                                           int *weights, int sumweights,
                                           int *subset, int Nsubset,    
                                           double *Q_ans);
void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights,
                                   int sumweights, double *ExpInf, 
                                   double *QQ_sym_ans) ;
void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights,
                                          int sumweights, int *subset,
                                          int Nsubset, double *ExpInf,
                                          double *QQ_sym_ans);
void C_VarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                 int sumweights, double *ExpInf, double *Q_ans);
void C_ExpectationX_weights(double* x, int N, int P, int *weights, double *P_ans);
void C_CovarianceX_weights(double *x, int N, int P, 
                           int *weights, double *PP_sym_ans);
void C_VarianceX_weights(double *x, int N, int P, int *weights, double *P_ans);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        