
#include "libcoin.h"
#include "helpers.h"

void C_LinearStatistic(double *x, int N, int P, double *y, int Q, double *PQ_ans) {
     C_KronSums(x, N, P, y, Q, PQ_ans);
}
     
void C_LinearStatistic_weights(double *x, int N, int P, double *y, int Q, 
                               double *weights, double *PQ_ans) {
     C_KronSums_weights(x, N, P, y, Q, weights, PQ_ans);
}
     
void C_LinearStatistic_subset(double *x, int N, int P, double *y, int Q, int *subset, 
                              int Nsubset, double *PQ_ans) {
     C_KronSums_subset(x, N, P, y, Q, subset, subset, Nsubset, PQ_ans);
}

void C_LinearStatistic_weights_subset(double *x, int N, int P, double *y, int Q, double *weights,
                                      int *subset, int Nsubset, double *PQ_ans) {
     C_KronSums_weights_subset(x, N, P, y, Q, weights, subset, Nsubset, PQ_ans);
}
     
void C_PermutedLinearStatistic(double *x, int N, int P, double *y, int Q, int *perm, 
                               int *original, int Nperm, double *PQ_ans) {
     C_KronSums_subset(x, N, P, y, Q, perm, original, Nperm, PQ_ans);
}

void C_ExpectationInfluence(double* y, int N, int Q, double *Q_ans) {
     C_colSums(y, N, Q, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / N;
}

void C_ExpectationInfluence_weights(double* y, int N, int Q, double *weights, 
                                    double sumweights, double *Q_ans) {
     C_colSums_weights(y, N, Q, weights, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}

void C_ExpectationInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                   double *Q_ans) {
     C_colSums_subset(y, N, Q, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / Nsubset;
}

void C_ExpectationInfluence_weights_subset(double* y, int N, int Q, double *weights, 
                                           double sumweights, int *subset, int Nsubset, 
                                           double *Q_ans) {
     C_colSums_weights_subset(y, N, Q, weights, subset, Nsubset, Q_ans);
     for (int q = 0; q < Q; q++) Q_ans[q] = Q_ans[q] / sumweights;
}
     
void C_CovarianceInfluence(double* y, int N, int Q, double *ExpInf, double *QQ_ans) {
     C_KronSums_center(y, N, Q, y, Q, ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / N;
}

void C_CovarianceInfluence_weights(double* y, int N, int Q, double *weights, 
                                   double sumweights, double *ExpInf, double *QQ_ans) {
     C_KronSums_weights_center(y, N, Q, y, Q, weights, ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / sumweights;
}

void C_CovarianceInfluence_subset(double* y, int N, int Q, int *subset, int Nsubset, 
                                  double *ExpInf, double *QQ_ans) {
     C_KronSums_subset_center(y, N, Q, y, Q, subset, subset, Nsubset, 
                              ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / Nsubset;
}

void C_CovarianceInfluence_weights_subset(double* y, int N, int Q, double *weights, 
                                          double sumweights, int *subset, int Nsubset, 
                                          double *ExpInf, double *QQ_ans) {
     C_KronSums_weights_subset_center(y, N, Q, y, Q, weights, subset, Nsubset, 
                                      ExpInf, ExpInf, QQ_ans);
     for (int q = 0; q < Q * Q; q++) QQ_ans[q] = QQ_ans[q] / sumweights;
}

void C_CovarianceX(double *x, int N, int P, double *PQ_ans) {
     C_KronSums(x, N, P, x, P, PQ_ans);
}
     
void C_CovarianceX_weights(double *x, int N, int P, 
                            double *weights, double *PQ_ans) {
     C_KronSums_weights(x, N, P, x, P, weights, PQ_ans);
}     
     
void C_CovarianceX_subset(double *x, int N, int P, int *subset, 
                          int Nsubset, double *PQ_ans) {
     C_KronSums_subset(x, N, P, x, P, subset, subset, Nsubset, PQ_ans);
}
     
void C_CovarianceX_weights_subset(double *x, int N, int P, double *weights,
                                  int *subset, int Nsubset, double *PQ_ans) {
     C_KronSums_weights_subset(x, N, P, x, P, weights, subset, Nsubset, PQ_ans);
}

void C_ExpectationLinearStatistic(int P, int Q, double *ExpInf, double *colSumsX, 
                                  double *PQ_ans) {

    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++)
            PQ_ans[q * P + p] = colSumsX[p] * ExpInf[q];
    }
}          

void C_CovarianceLinearStatistic(int P, int Q, double *CovInf, double *colSumsX, double *CovX, 
                                 double sumweights, double *PQPQ_tmp, double *PQQ_tmp, 
                                 double *PQPQ_ans) {

    double f1, f2;
    int PQ = P * Q;

    f1 = sumweights / (sumweights - 1);
    f2 = 1 / (sumweights - 1);
        
    if (PQ == 1) {
        PQPQ_ans[0] = f1 * CovInf[0] * CovX[0];
        PQPQ_ans[0] -= f2 * CovInf[0] * colSumsX[0] * colSumsX[0];
    } else {
        C_kronecker(CovInf, Q, Q, CovX, P, P, PQPQ_ans);
        C_kronecker(CovInf, Q, Q, colSumsX, P, 1, PQQ_tmp);
        C_kronecker(PQQ_tmp, PQ, Q, colSumsX, 1, P, PQPQ_tmp);  
 
        for (int k = 0; k < (PQ * PQ); k++)
             PQPQ_ans[k] = f1 * PQPQ_ans[k] - f2 * PQPQ_tmp[k];
    }
}
