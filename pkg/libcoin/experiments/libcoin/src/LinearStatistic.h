
/* Variables
   x:		a double N x P matrix
   y:		a double N x Q matrix / a double M x Q matrix
   weights:	an integer N vector with sumweights = sum(weights)
   subset:	an integer Nsubset vector
   ExpInf:	expectation of influence function y, a double Q vector
   CovInf:	covariance of influence function y, a double Q x Q matrix
   CovInf:	variance of influence function y, a double Q vector
   ExpX:	"expectation" of x, a double P vector
   CovX:	"covariance" of x, a double P x P matrix
   VarX:	"variance" of x, a double P vector
   PQ_ans:	return value, a double P x Q matrix
   PP_ans:	return value, a double P x Q matrix
   QQ_ans:	return value, a double Q x Q matrix
   P_ans:	return value, a double P vector
   Q_ans:	return value, a double Q vector
   PP_tmp:	temp variable, a double P x P matrix
   P_tmp:	temp variable, a double P vector
*/   
   
void C_LinearStatistic(double *x, int N, int P, double *y, int Q, 
                       double *PQ_ans); 
void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans);
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, 
                              int *subset, int Nsubset, double *PQ_ans); 
void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, 
                                      int *weights, int *subset, int Nsubset, 
                                      double *PQ_ans); 
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, 
                               int *perm, int *original, int Nperm, 
                               double *PQ_ans); 
void C_LinearStatistic_2d(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans); 
void C_LinearStatistic_maxstat(int *ix, int N, int P, double *y, int Q, 
                               double *PQ_ans);
void C_LinearStatistic_maxstat_weights(int *ix, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans);
void C_LinearStatistic_maxstat_subset(int *ix, int N, int P, double *y, int Q, 
                              int *subset, int Nsubset, double *PQ_ans) ;
void C_LinearStatistic_maxstat_weights_subset(int *ix, int N, int P, double *y, int Q, 
                                      int *weights, int *subset, int Nsubset, 
                                      double *PQ_ans);
void C_LinearStatistic_maxstat_2d(int N, int P, double *y, int M, int Q, 
                                  int *weights, double *PQ_ans);
void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans); 
void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights, 
                                    int sumweights, double *Q_ans); 
void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, 
                                   int Nsubset, double *Q_ans); 
void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, 
                                           int *weights, int sumweights, 
                                           int *subset, int Nsubset, 
                                           double *Q_ans); 
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, 
                           double *QQ_ans); 
void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, 
                                   double *QQ_ans); 
void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                  int Nsubset, double *ExpInf, double *QQ_ans); 
void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, 
                                          int Nsubset, double *ExpInf, 
                                          double *QQ_ans); 
void C_VarianceInfluence(double* y, int N, int Q, double *ExpInf, double *Q_ans); 
void C_VarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                 int sumweights, double *ExpInf, double *Q_ans);
void C_VarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                int Nsubset, double *ExpInf, double *Q_ans); 
void C_VarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                        int sumweights, int *subset, 
                                        int Nsubset, double *ExpInf, 
                                        double *Q_ans); 
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
void C_VarianceX(double *x, int N, int P, double *P_ans); 
void C_VarianceX_weights(double *x, int N, int P, int *weights, double *P_ans); 
void C_VarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *P_ans); 
void C_VarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *P_ans); 
void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *ExpX, 
                                  double *PQ_ans);
void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *ExpX, 
                                 double *CovX, int sumweights, double *PP_tmp, int add,
                                 double *PQPQ_ans); 
void C_VarianceLinearStatistic(int P, int Q, double *VarInf, double *ExpX, 
                               double *VarX, int sumweights, double *P_tmp, int add,
                               double *PQ_ans); 
