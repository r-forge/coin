
#include "Sums.h"
#include "helpers.h"

void C_LinearStatistic(double *x, int N, int P, double *y, int Q, double *PQ_ans) {
     C_KronSums(x, N, P, y, Q, PQ_ans);
}
     
void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *PQ_ans) {
     C_KronSums_weights(x, N, P, y, Q, weights, PQ_ans);
}
     
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, int *subset, 
                              int Nsubset, double *PQ_ans) {
     C_KronSums_subset(x, N, P, y, Q, subset, subset, Nsubset, PQ_ans);
}

void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, int *weights,
                                      int *subset, int Nsubset, double *PQ_ans) {
     C_KronSums_weights_subset(x, N, P, y, Q, weights, subset, Nsubset, PQ_ans);
}
     
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, int *perm, 
                               int *original, int Nperm, double *PQ_ans) {
     C_KronSums_subset(x, N, P, y, Q, perm, original, Nperm, PQ_ans);
}

void C_LinearStatistic_2d(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans) {
    C_KronSums_2dweights(x, N, P, y, M, Q, weights, PQ_ans);
}

void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans) {
     C_colSums(y, N, Q, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_ExpectationInfluence_weights(double* y, int N, int Q, int *weights, 
                                    int sumweights, double *Q_ans) {
     C_colSums_weights(y, N, Q, weights, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                   double *Q_ans) {
     C_colSums_subset(y, N, Q, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                           int sumweights, int *subset, int Nsubset, 
                                           double *Q_ans) {
     C_colSums_weights_subset(y, N, Q, weights, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}
     
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, double *QQ_ans) {
     C_KronSums_center(y, N, Q, y, Q, ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / N;
}

void C_CovarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, double *QQ_ans) {
     C_KronSums_weights_center(y, N, Q, y, Q, weights, ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / sumweights;
}

void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                  double *ExpInf, double *QQ_ans) {
     C_KronSums_subset_center(y, N, Q, y, Q, subset, subset, Nsubset, 
                              ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / Nsubset;
}

void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, int Nsubset, 
                                          double *ExpInf, double *QQ_ans) {
     C_KronSums_weights_subset_center(y, N, Q, y, Q, weights, subset, Nsubset, 
                                      ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / sumweights;
}

void C_VarianceInfluence(double* y, int N, int Q, double *ExpInf, double *Q_ans) {
     C_colSums2_center(y, N, Q, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_VarianceInfluence_weights(double* y, int N, int Q, int *weights, 
                                   int sumweights, double *ExpInf, double *Q_ans) {
     C_colSums2_weights_center(y, N, Q, weights, ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_VarianceInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                  double *ExpInf, double *Q_ans) {
     C_colSums2_subset_center(y, N, Q, subset, Nsubset, 
                              ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_VarianceInfluence_weights_subset(double* y, int N, int Q, int *weights, 
                                          int sumweights, int *subset, int Nsubset, 
                                          double *ExpInf, double *Q_ans) {
     C_colSums2_weights_subset_center(y, N, Q, weights, subset, Nsubset, 
                                      ExpInf, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationX(double* x, int N, int P, double *P_ans) {
     C_colSums(x, N, P, P_ans);
}

void C_ExpectationX_weights(double* x, int N, int P, int *weights, 
                            int sumweights, double *P_ans) {
     C_colSums_weights(x, N, P, weights, P_ans);
}

void C_ExpectationX_subset(double* x, int N, int P, int *subset, int Nsubset, 
                           double *P_ans) {
     C_colSums_subset(x, N, P, subset, Nsubset, P_ans);
}

void C_ExpectationX_weights_subset(double* x, int N, int P, int *weights, 
                                   int sumweights, int *subset, int Nsubset, 
                                   double *P_ans) {
     C_colSums_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_CovarianceX(double *x, int N, int P, double *PP_ans) {
     C_KronSums(x, N, P, x, P, PP_ans);
}
     
void C_CovarianceX_weights(double *x, int N, int P, 
                            int *weights, double *PP_ans) {
     C_KronSums_weights(x, N, P, x, P, weights, PP_ans);
}     
     
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PP_ans) {
     C_KronSums_subset(x, N, P, x, P, subset, subset, Nsubset, PP_ans);
}
     
void C_CovarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *PP_ans) {
     C_KronSums_weights_subset(x, N, P, x, P, weights, subset, Nsubset, PP_ans);
}

void C_VarianceX(double *x, int N, int P, double *P_ans) {
     C_colSums2(x, N, P, P_ans);
}
     
void C_VarianceX_weights(double *x, int N, int P, 
                            int *weights, double *P_ans) {
     C_colSums2_weights(x, N, P, weights, P_ans);
}     
     
void C_VarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *P_ans) {
     C_colSums2_subset(x, N, P, subset, Nsubset, P_ans);
}
     
void C_VarianceX_weights_subset(double *x, int N, int P, int *weights,
                                  int *subset, int Nsubset, double *P_ans) {
     C_colSums2_weights_subset(x, N, P, weights, subset, Nsubset, P_ans);
}

void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *ExpX, 
                                  double *PQ_ans) {

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] = ExpX[p] * ExpInf[q];
    }
}          

void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *ExpX, double *CovX, 
                                 int sumweights, double *PP_tmp, double *PQPQ_ans) {

    double f1, f2;
    int PQ = P * Q;

    f1 = (double) sumweights / (sumweights - 1);
    f2 = 1.0 / (sumweights - 1);
        
    if (PQ == 1) {
        PQPQ_ans[0] = f1 * CovInf[0] * CovX[0];
        PQPQ_ans[0] -= f2 * CovInf[0] * ExpX[0] * ExpX[0];
    } else {
        C_kronecker(ExpX, P, 1, ExpX, 1, P, PP_tmp);
        for (int p = 0; p < P * P; p++)
            PP_tmp[p] = f1 * CovX[p] - f2 * PP_tmp[p];
        C_kronecker(CovInf, Q, Q, PP_tmp, P, P, PQPQ_ans);
    }
}

void C_VarianceLinearStatistic(int P, int Q, double *VarInf, double *ExpX, double *VarX, 
                               int sumweights, double *P_tmp, double *PQ_ans) {

    double f1, f2;
    int PQ = P * Q;

    if (PQ == 1) {
        C_CovarianceLinearStatistic(P, Q, VarInf, ExpX, VarX, 
                                    sumweights, P_tmp, PQ_ans);
    } else {

        f1 = (double) sumweights / (sumweights - 1);
        f2 = 1.0 / (sumweights - 1);
        
        for (int p = 0; p < P; p++)
            P_tmp[p] = f1 * VarX[p] - f2 * ExpX[p] * ExpX[p];
        
        C_kronecker(VarInf, 1, Q, P_tmp, 1, P, PQ_ans);
    }
}

