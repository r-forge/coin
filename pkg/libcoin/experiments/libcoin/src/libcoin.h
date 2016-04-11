
#include <R.h>
#include <Rinternals.h>
/*
void C_LinearStatistic(double *x, int N, int P, double *y, int Q, double *PQ_ans);
void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               double *weights, double *PQ_ans) { 
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, int *subset, 
                              int Nsubset, double *PQ_ans);
void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, double *weights,
                                      int *subset, int Nsubset, double *PQ_ans);
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, int *perm, 
                               int *original, int Nperm, double *PQ_ans);
void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans);
void C_ExpectationInfluence_weights(double* y, int N, int Q, double *weights, 
                                    double sumweights, double *Q_ans);
void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                   double *Q_ans);
void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, double *weights, 
                                           double sumweights, int *subset, int Nsubset, 
                                           double *Q_ans);
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_weights(double* y, int N, int Q, double *weights, 
                                   double sumweights, double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                  double *ExpInf, double *QQ_ans);
void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, double *weights, 
                                          double sumweights, int *subset, int Nsubset, 
                                          double *ExpInf, double *QQ_ans);
void C_CovarianceX(double *x, int N, int P, double *PQ_ans);
void C_CovarianceX_weights(double *x, int N, int P, 
                            double *weights, double *PQ_ans);
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PQ_ans);
void C_CovarianceX_weights_subset(double *x, int N, int P, double *weights,
                                  int *subset, int Nsubset, double *PQ_ans);
void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *colSumsX, 
                                  double *PQ_ans);
void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *colSumsX, double *CovX, 
                                 double sumweights, double *PQPQ_tmp, double *PQQ_tmp, 
                                 double *PQPQ_ans);
*/
void C_sum(double *x, int N, double *ans);
void C_sum_subset(double *x, int N, int *subset, int Nsubset, double *ans);
void C_colSums(double *x, int N, int P, double *ans);
void C_colSums_weights(double *x, int N, int P, double *weights, double *ans);
void C_colSums_subset(double *x, int N, int P, int *subset, int Nsubset, double *ans);
void C_colSums_weights_subset(double *x, int N, int P, double *weights, 
                              int *subset, int Nsubset, double *ans);
void C_KronSums(double *x, int N, int P, double *y, int Q, double *ans);
void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        double *weights, double *ans);
void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, double *ans);
void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               double *weights, int *subset, int Nsubset, 
                               double *ans);
void C_KronSums_2dweights(double *x, int N, int P, double *y, int M, int Q, 
                          double *weights, double *ans);
void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *ans);
void C_KronSums_weights_center(double *x, int N, int P, double *y, int Q, 
                               double *weights, double *centerx, 
                               double *centery, double *ans);
void C_KronSums_subset_center(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *ans);
void C_KronSums_weights_subset_center(double *x, int N, int P, double *y, 
                                      int Q, double *weights, 
                                      int *subset, int Nsubset, 
                                      double *centerx, double *centery, double *ans);
