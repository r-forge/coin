
#include <R.h>
#include <Rinternals.h>

/* sum(x) */

void C_sum(double *x, int N, double *ans) {

    ans[0] = 0.0;
    for (int i; i < N; i++) ans[0] += x[i];
}

/* sum(x[subset]) */

void C_sum_subset(double *x, int N, int *subset, int Nsubset, double *ans) {

    ans[0] = 0.0;
    for (int i; i < Nsubset; i++) ans[0] += x[subset[i]];
}

/* colSums(x) */

void C_colSums(double *x, int N, int P, double *ans) {

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int i = 0; i < N; i++) {
        for (int p = 0; p < P; p++)
            ans[p] += x[p * N + i];
    }
}

/* colSums(x * w) */

void C_colSums_weights(double *x, int N, int P, int *weights, double *ans) {

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int i = 0; i < N; i++) {
        for (int p = 0; p < P; p++)
            ans[p] += weights[i] * x[p * N + i];
    }
}

/* colSums(x[subset,]) */

void C_colSums_subset(double *x, int N, int P, int *subset, int Nsubset, double *ans) {

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int i = 0; i < Nsubset; i++) {
        for (int p = 0; p < P; p++)
            ans[p] += x[p * N + subset[i]];
    }
}

/* colSums(x[subset,] * weights[subset]) */

void C_colSums_weights_subset(double *x, int N, int P, int *weights, 
                              int *subset, int Nsubset, double *ans) {

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int i = 0; i < Nsubset; i++) {
        for (int p = 0; p < P; p++)
            ans[p] += weights[subset[i]] * x[p * N + subset[i]];
    }
}

/* sum_i (t(x[i,]) %*% y[i,]) */

void C_KronSums(double *x, int N, int P, double *y, int Q, double *ans) {

    int qP, qN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             for (int p = 0; p < P; p++)
                 ans[qP + p] +=  y[qN + i] * x[p * N + i];
        }
    }
}

/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */

void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        int *weights, double *ans) {

    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = y[qN + i] * weights[i];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += x[p * N + i] * tmp;
        }
    }
}

void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, double *ans) {

    int qP, qN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             for (int p = 0; p < P; p++)
                 ans[qP + p] += y[qN + subsety[i]] * x[p * N + subsetx[i]];
        }
    }
}

void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               int *weights, int *subset, int Nsubset, 
                               double *ans) {
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = y[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += x[p * N + subset[i]] * tmp;
        }
    }
}

void C_KronSums_2dweights(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *ans) {

    int qPp, qM, pNi, iM;
        
    for (int p = 0; p < P * Q; p++) ans[p] = 0.0;
        
    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++) {
            qPp = q * P + p;
            qM = q * M;
            for (int i = 0; i < N; i++) {
                pNi = p * N + i;
                iM = i * M;
                for (int m = 0; m < M; m++)
                      ans[qPp] += y[qM + m] * x[pNi] * weights[iM + m];
            }
        }
    }
}

/* sum_i (t(x[i,] - centerx) %*% (y[i,] - centery)) */

void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *ans) {

    int qP, qN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             for (int p = 0; p < P; p++)
                 ans[qP + p] +=  y[qN + i] * x[p * N + i];
        }
    }
}

void C_KronSums_weights_center(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *centerx, 
                               double *centery, double *ans) {

    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = y[qN + i] * weights[i];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += x[p * N + i] * tmp;
        }
    }
}

void C_KronSums_subset_center(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *ans) {

    int qP, qN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             for (int p = 0; p < P; p++)
                 ans[qP + p] += y[qN + subsety[i]] * x[p * N + subsetx[i]];
        }
    }
}

void C_KronSums_weights_subset_center(double *x, int N, int P, double *y, 
                                      int Q, int *weights, 
                                      int *subset, int Nsubset, 
                                      double *centerx, double *centery, double *ans) {
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = y[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += x[p * N + subset[i]] * tmp;
        }
    }
}
