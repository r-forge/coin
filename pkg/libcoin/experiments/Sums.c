
#include <R.h>
#include <Rinternals.h>

void C_sum(double *x, int n, double *ans) {

    ans[0] = 0.0;
    for (int i; i < n; i++) ans[0] += x[i];
}

void C_sum_subset(double *x, int n, int *subset, int nsubset, double *ans) {

    ans[0] = 0.0;
    for (int i; i < nsubset; i++) ans[0] += x[subset[i]];
}

void C_colSums(double *x, int n, int p, double *ans) {

    for (int j = 0; j < p; j++)
        ans[j] = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++)
            ans[j] += x[j * n + i];
    }
}

void C_colSums_weights(double *x, int n, int p, double *w, double *ans) {

    for (int j = 0; j < p; j++)
        ans[j] = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++)
            ans[j] += w[i] * x[j * n + i];
    }
}

void C_colSums_subset(double *x, int n, int p, int *subset, int nsubset, double *ans) {

    for (int j = 0; j < p; j++)
        ans[j] = 0.0;

    for (int i = 0; i < nsubset; i++) {
        for (int j = 0; j < p; j++)
            ans[j] += x[j * n + subset[i]];
    }
}
void C_colSums_weights_subset(double *x, int n, int p, double *w, 
                              int *subset, int nsubset, double *ans) {

    for (int j = 0; j < p; j++)
        ans[j] = 0.0;

    for (int i = 0; i < nsubset; i++) {
        for (int j = 0; j < p; j++)
            ans[j] += w[subset[i]] * x[j * n + subset[i]];
    }
}

void C_KronSums(double *x, int n, int p, double *y, int q, double *ans) {

    int kp, kn;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < n; i++) {
             for (int j = 0; j < p; j++)
                 ans[kp + j] +=  y[kn + i] * x[j*n + i];
        }
    }
}


void C_KronSums_weights(double *x, int n, int p, double *y, int q, double *w, double *ans) {

    int kp, kn;
    double tmp;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < n; i++) {
             tmp = y[kn + i] * w[i];
             for (int j = 0; j < p; j++)
                 ans[kp + j] += x[j*n + i] * tmp;
        }
    }
}

void C_KronSums_subset(double *x, int n, int p, double *y, int q, 
                       int *subset, int nsubset, double *ans) {

    int kp, kn;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < nsubset; i++) {
             for (int j = 0; j < p; j++)
                 ans[kp + j] += y[kn + subset[i]] * x[j*n + subset[i]];
        }
    }
}

void C_KronSums_weights_subset(double *x, int n, int p, double *y, int q, double *w, 
                               int *subset, int nsubset, double *ans) {
    int kp, kn;
    double tmp;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < nsubset; i++) {
             tmp = y[kn + subset[i]] * w[subset[i]];
             for (int j = 0; j < p; j++)
                 ans[kp + j] += x[j*n + subset[i]] * tmp;
        }
    }
}

void C_KronSums_center(double *x, int n, int p, double *y, int q, 
                       double *centerx, double *centery, double *ans) {

    int kp, kn;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < n; i++) {
             for (int j = 0; j < p; j++)
                 ans[kp + j] +=  y[kn + i] * x[j*n + i];
        }
    }
}

void C_KronSums_weights_center(double *x, int n, int p, double *y, int q, double *w, 
                               double *centerx, double *centery, double *ans) {

    int kp, kn;
    double tmp;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < n; i++) {
             tmp = y[kn + i] * w[i];
             for (int j = 0; j < p; j++)
                 ans[kp + j] += x[j*n + i] * tmp;
        }
    }
}

void C_KronSums_subset_center(double *x, int n, int p, double *y, int q, int *subset, 
                              int nsubset, double *centerx, double *centery, double *ans) {

    int kp, kn;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < nsubset; i++) {
             for (int j = 0; j < p; j++)
                 ans[kp + j] += y[kn + subset[i]] * x[j*n + subset[i]];
        }
    }
}

void C_KronSums_weights_subset_center(double *x, int n, int p, double *y, int q, double *w, 
                                      int *subset, int nsubset, 
                                      double *centerx, double *centery, double *ans) {
    int kp, kn;
    double tmp;
        
    for (int k = 0; k < q; k++) {
        kn = k * n;
        kp = k * p;
        for (int j = 0; j < p; j++) ans[kp + j] = 0.0;
        for (int i = 0; i < nsubset; i++) {
             tmp = y[kn + subset[i]] * w[subset[i]];
             for (int j = 0; j < p; j++)
                 ans[kp + j] += x[j*n + subset[i]] * tmp;
        }
    }
}
