
void C_LinearStatistic(double *x, int N, int P, double *y, int Q, double *PQ_ans);
void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans);
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, int *subset, 
                              int Nsubset, double *PQ_ans);
void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, int *weights,
                                      int *subset, int Nsubset, double *PQ_ans);
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, int *perm, 
                               int *original, int Nperm, double *PQ_ans);
void C_LinearStatistic_2d(double *x, int N, int P, double *y, int M, int Q,
                          int *weights, double *PQ_ans);
void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans);
void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights, 
                                    int sumweights, double *Q_ans);
void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                   double *Q_ans);
void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                           int sumweights, int *subset, int Nsubset, 
                                           double *Q_ans);
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                  double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, int Nsubset, 
                                          double *ExpInf, double *QQ_ans);
void C_ExpectationX(double* x, int N, int P, double *P_ans);
void C_ExpectationX_weights(double* x, int N, int P, int *weights, 
                                    int sumweights, double *P_ans);
void C_ExpectationX_subset(double* x, int N, int P, int *subset, int Nsubset, 
                                   double *P_ans);
void C_ExpectationX_weights_subset(double* x, int N, int P, int *weights, 
                                           int sumweights, int *subset, int Nsubset, 
                                           double *P_ans);
void C_CovarianceX(double *x, int N, int P, double *PP_ans);
void C_CovarianceX_weights(double *x, int N, int P, 
                            int *weights, double *PP_ans);
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PP_ans);
void C_CovarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *PP_ans);
void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *colSumsX, 
                                  double *PQ_ans);
void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *colSumsX, double *CovX, 
                                 int sumweights, double *PQPQ_tmp, double *PQQ_tmp, 
                                 double *PQPQ_ans);
