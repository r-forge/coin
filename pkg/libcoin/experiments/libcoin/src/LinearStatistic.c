
#include "Sums.h"
#include "helpers.h"


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
                       double *PQ_ans) 
{
     C_KronSums(x, N, P, y, Q, PQ_ans);
}

void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans)
{                               
     C_KronSums_weights(x, N, P, y, Q, weights, PQ_ans);
}
     
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, 
                              int *subset, int Nsubset, double *PQ_ans) 
{
     C_KronSums_subset(x, N, P, y, Q, subset, subset, Nsubset, PQ_ans);
}

void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, 
                                      int *weights, int *subset, int Nsubset, 
                                      double *PQ_ans) 
{
     C_KronSums_weights_subset(x, N, P, y, Q, weights, subset, Nsubset, PQ_ans);
}
     
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, 
                               int *perm, int *original, int Nperm, 
                               double *PQ_ans) 
{
     C_KronSums_subset(x, N, P, y, Q, perm, original, Nperm, PQ_ans);
}

void C_LinearStatistic_2d(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans) 
{
    C_KronSums_2dweights(x, N, P, y, M, Q, weights, PQ_ans);
}

void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans) 
{
     C_colSums(y, N, Q, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights, 
                                    int sumweights, double *Q_ans) 
{
     C_colSums_weights(y, N, Q, weights, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, 
                                   int Nsubset, double *Q_ans) 
{
     C_colSums_subset(y, N, Q, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, 
                                           int *weights, int sumweights, 
                                           int *subset, int Nsubset, 
                                           double *Q_ans) 
{
     C_colSums_weights_subset(y, N, Q, weights, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}
     
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, 
                           double *QQ_sym_ans) 
{
     C_KronSums_sym_center(y, N, Q, ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) QQ_sym_ans[q] = QQ_sym_ans[q] / N;
}

void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, 
                                   double *QQ_sym_ans) 
{
     C_KronSums_sym_center_weights(y, N, Q, weights, ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
}

void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                  int Nsubset, double *ExpInf, double *QQ_sym_ans) 
{
     C_KronSums_sym_center_subset(y, N, Q, subset, Nsubset, 
                              ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) QQ_sym_ans[q] = QQ_sym_ans[q] / Nsubset;
}

void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, 
                                          int Nsubset, double *ExpInf, 
                                          double *QQ_sym_ans) 
{
     C_KronSums_sym_center_weights_subset(y, N, Q, weights, subset, Nsubset, 
                                      ExpInf, QQ_sym_ans);
     for (int q = 0; q < Q * (Q + 1) / 2; q++) QQ_sym_ans[q] = QQ_sym_ans[q] / sumweights;
}

void C_VarianceInfluence(double* y, int N, int Q, double *ExpInf, double *Q_ans) 
{
     C_colSums2_center(y, N, Q, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_VarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                 int sumweights, double *ExpInf, double *Q_ans)
{
     C_colSums2_center_weights(y, N, Q, weights, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_VarianceInfluence_subset(double* y, int N, int Q, int *subset, 
                                int Nsubset, double *ExpInf, double *Q_ans) 
{
     C_colSums2_center_subset(y, N, Q, subset, Nsubset, 
                              ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_VarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                        int sumweights, int *subset, 
                                        int Nsubset, double *ExpInf, 
                                        double *Q_ans) 
{
     C_colSums2_center_weights_subset(y, N, Q, weights, subset, Nsubset, 
                                      ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationX(double* x, int N, int P, double *P_ans) 
{
     C_colSums(x, N, P, P_ans);
}

void C_ExpectationX_weights(double* x, int N, int P, int *weights, 
                            int sumweights, double *P_ans) 
{
     C_colSums_weights(x, N, P, weights, P_ans);
}

void C_ExpectationX_subset(double* x, int N, int P, int *subset, int Nsubset, 
                           double *P_ans) 
{
     C_colSums_subset(x, N, P, subset, Nsubset, P_ans);
}

void C_ExpectationX_weights_subset(double* x, int N, int P, int *weights, 
                                   int sumweights, int *subset, int Nsubset, 
                                   double *P_ans) 
{
     C_colSums_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_CovarianceX(double *x, int N, int P, double *PP_sym_ans) 
{
     C_KronSums_sym(x, N, P, PP_sym_ans);
}
     
void C_CovarianceX_weights(double *x, int N, int P, 
                           int *weights, double *PP_sym_ans) 
{
     C_KronSums_sym_weights(x, N, P, weights, PP_sym_ans);
}     
     
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PP_sym_ans) 
{
     C_KronSums_sym_subset(x, N, P, subset, Nsubset, PP_sym_ans);
}
     
void C_CovarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *PP_sym_ans) 
{
     C_KronSums_sym_weights_subset(x, N, P, weights, subset, Nsubset, PP_sym_ans);
}

void C_VarianceX(double *x, int N, int P, double *P_ans) 
{
     C_colSums2(x, N, P, P_ans);
}
     
void C_VarianceX_weights(double *x, int N, int P, int *weights, double *P_ans) 
{
     C_colSums2_weights(x, N, P, weights, P_ans);
}     
     
void C_VarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *P_ans) 
{
     C_colSums2_subset(x, N, P, subset, Nsubset, P_ans);
}
     
void C_VarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *P_ans) 
{
     C_colSums2_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *ExpX, 
                                  double *PQ_ans)
{
    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] = ExpX[p] * ExpInf[q];
    }
}          

void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *ExpX, 
                                 double *CovX, int sumweights, double *PP_sym_tmp, 
                                 double *PQPQ_sym_ans) 
{
    double f1 = (double) sumweights / (sumweights - 1);
    double f2 = 1.0 / (sumweights - 1);
        
    if (P * Q == 1) {
        PQPQ_sym_ans[0] = f1 * CovInf[0] * CovX[0];
        PQPQ_sym_ans[0] -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
    } else {
        C_KronSums_sym(ExpX, 1, P, PP_sym_tmp);
        for (int p = 0; p < P * (P + 1) / 2; p++)
            PP_sym_tmp[p] = f1 * CovX[p] - f2 * PP_sym_tmp[p];
        C_kronecker_sym(CovInf, Q, PP_sym_tmp, P, PQPQ_sym_ans);
    }
}

void C_VarianceLinearStatistic(int P, int Q, double *VarInf, double *ExpX, 
                               double *VarX, int sumweights, double *P_tmp, 
                               double *PQ_ans) 
{
    if (P * Q == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX, 
                                    sumweights, P_tmp, PQ_ans);
    } else {

        double f1 = (double) sumweights / (sumweights - 1);
        double f2 = 1.0 / (sumweights - 1);
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, PQ_ans);
    }
}
