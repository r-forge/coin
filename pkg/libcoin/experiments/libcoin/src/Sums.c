
#include <R.h>
#include <Rinternals.h>

/* sum(x) */

int C_sum(int *x, int N) {

    int ans = 0;
    for (int i; i < N; i++) ans += x[i];
    return(ans);
}

/* sum(x[subset]) */

int C_sum_subset(int *x, int N, int *subset, int Nsubset) {

    int ans = 0;
    for (int i; i < Nsubset; i++) ans += x[subset[i]];
    return(ans);
}

/* colSums(x) */

void C_colSums(double *x, int N, int P, double *ans) {

    int pN;

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        pN = p * N;
        for (int i = 0; i < N; i++)
            ans[p] += x[pN + i];
    }
}

/* colSums(x * w) */

void C_colSums_weights(double *x, int N, int P, int *weights, double *ans) {

    int pN;

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        pN = p * N;
        for (int i = 0; i < N; i++)
            ans[p] += weights[i] * x[pN + i];
    }
}

/* colSums(x[subset,]) */

void C_colSums_subset(double *x, int N, int P, int *subset, int Nsubset, double *ans) {

    int pN;

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            ans[p] += x[pN + subset[i]];
    }
}

/* colSums(x[subset,] * weights[subset]) */

void C_colSums_weights_subset(double *x, int N, int P, int *weights, 
                              int *subset, int Nsubset, double *ans) {

    int pN;

    for (int p = 0; p < P; p++)
        ans[p] = 0.0;

    for (int p = 0; p < P; p++) {
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            ans[p] += weights[subset[i]] * x[pN + subset[i]];
    }
}

/* sum_i (t(x[i,]) %*% y[i,]) */

void C_KronSums(double *x, int N, int P, double *y, int Q, double *ans) {

    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int p = 0; p < P; p++) {
            pN = p * N;
            for (int i = 0; i < N; i++)
                 ans[qP + p] +=  y[qN + i] * x[pN + i];
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

    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int p = 0; p < P; p++) {
             pN = p * N;
        for (int i = 0; i < Nsubset; i++)
                 ans[qP + p] += y[qN + subsety[i]] * x[pN + subsetx[i]];
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

    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int p = 0; p < P; p++) {
            pN = p * N;
            for (int i = 0; i < N; i++)
                 ans[qP + p] +=  (y[qN + i] - centery[q]) * 
                                 (x[pN + i] - centerx[p]);
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
             tmp = (y[qN + i] - centery[q]) * weights[i];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += (x[p * N + i] - centerx[p]) * tmp;
        }
    }
}

void C_KronSums_subset_center(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *ans) {

    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) ans[qP + p] = 0.0;
        for (int p = 0; p < P; p++) {
            pN = p * N;
            for (int i = 0; i < Nsubset; i++)
                 ans[qP + p] += (y[qN + subsety[i]] - centery[q]) * 
                                (x[pN + subsetx[i]] - centerx[p]);
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
             tmp = (y[qN + subset[i]] - centery[q]) * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 ans[qP + p] += (x[p * N + subset[i]] - centerx[p]) * tmp;
        }
    }
}
