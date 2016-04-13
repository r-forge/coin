
#include <R.h>
#include <Rinternals.h>

/* Variables
   x:           a double N x P matrix
   y:           a double N x Q matrix / a double M x Q matrix
   weights:     an integer N vector with sumweights = sum(weights)
   subsetx:     an integer Nsubset vector
   subsety:     an integer Nsubset vector
   centerx:	a double P vector centering the columns of x
   centery:	a double Q vector centering the columns of y
   PQ_ans:      return value, a double P x Q matrix
   P_ans:       return value, a double P vector
   Q_ans:       return value, a double Q vector
*/
                                                   

/* sum(weights) */
int C_sum(int *weights, int N)
{
    int ans = 0;
    for (int i; i < N; i++) ans += weights[i];
    return(ans);
}

/* sum(weights[subset]) */
int C_sum_subset(int *weights, int N, int *subset, int Nsubset)
{
    int ans = 0;
    for (int i; i < Nsubset; i++) ans += weights[subset[i]];
    return(ans);
}

/* colSums(x) */
void C_colSums(double *x, int N, int P, double *P_ans) 
{
    int pN;
    
    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += x[pN + i];
    }
}

/* colSums(x * w) */
void C_colSums_weights(double *x, int N, int P, int *weights, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * x[pN + i];
    }
}

/* colSums(x[subsetx,]) */
void C_colSums_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                      double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += x[pN + subsetx[i]];
    }
}

/* colSums(x[subsetx,] * weights[subsetx]) */
void C_colSums_weights_subset(double *x, int N, int P, int *weights, 
                              int *subsetx, int Nsubset, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] * x[pN + subsetx[i]];
    }
}

/* colSums(x^2) */
void C_colSums2(double *x, int N, int P, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += pow(x[pN + i], 2);
    }
}

/* colSums(x^2 * w) */
void C_colSums2_weights(double *x, int N, int P, int *weights, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * pow(x[pN + i], 2);
    }
}

/* colSums(x[subsetx,]^2) */
void C_colSums2_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                       double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += pow(x[pN + subsetx[i]], 2);
    }
}

/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_weights_subset(double *x, int N, int P, int *weights, 
                               int *subsetx, int Nsubset, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] * pow(x[pN + subsetx[i]], 2);
    }
}


/* colSums((x-center)^2) */
void C_colSums2_center(double *x, int N, int P, double *centerx, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += pow(x[pN + i] - centerx[p], 2);
    }
}

/* colSums((x-center)^2 * w) */
void C_colSums2_center_weights(double *x, int N, int P, int *weights, 
                               double *centerx, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += weights[i] * pow(x[pN + i] - centerx[p], 2);
    }
}

/* colSums((x[subsetx,] - center)^2) */
void C_colSums2_center_subset(double *x, int N, int P, int *subsetx, 
                              int Nsubset, double *centerx, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += pow(x[pN + subsetx[i]] - centerx[p], 2);
    }
}

/* colSums((x[subsetx,]-center)^2 * weights[subsetx]) */
void C_colSums2_center_weights_subset(double *x, int N, int P, int *weights, 
                                      int *subsetx, int Nsubset, 
                                      double *centerx, double *P_ans) 
{
    int pN;

    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.0;
        pN = p * N;
        for (int i = 0; i < Nsubset; i++)
            P_ans[p] += weights[subsetx[i]] * 
                        pow(x[pN + subsetx[i]] - centerx[p], 2);
    }
}


/* sum_i (t(x[i,]) %*% y[i,]) */
void C_KronSums(double *x, int N, int P, double *y, int Q, double *PQ_ans) 
{
    int pN, qP, qN;
    
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            for (int i = 0; i < N; i++)
                 PQ_ans[qP + p] +=  y[qN + i] * x[pN + i];
        }
    }
}

/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        int *weights, double *PQ_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = y[qN + i] * weights[i];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += x[p * N + i] * tmp;
        }
    }
}

void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, double *PQ_ans) 
{
    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            for (int i = 0; i < Nsubset; i++)
                PQ_ans[qP + p] += y[qN + subsety[i]] * x[pN + subsetx[i]];
        }
    }
}

void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               int *weights, int *subset, int Nsubset, 
                               double *PQ_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = y[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += x[p * N + subset[i]] * tmp;
        }
    }
}

void C_KronSums_2dweights(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans) 
{
    int qPp, qM, pNi;
        
    for (int p = 0; p < P; p++) {
        for (int q = 0; q < Q; q++) {
            PQ_ans[q * P + p] = 0.0;
            qPp = q * P + p;
            qM = q * M;
            for (int i = 0; i < N; i++) {
                pNi = p * N + i;
                for (int m = 0; m < M; m++)
                      PQ_ans[qPp] += y[qM + m] * x[pNi] * weights[m * N + i];
            }
        }
    }
}

/* sum_i (t(x[i,] - centerx) %*% (y[i,] - centery)) */
void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *PQ_ans) 
{
    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            for (int i = 0; i < N; i++)
                 PQ_ans[qP + p] +=  (y[qN + i] - centery[q]) * 
                                    (x[pN + i] - centerx[p]);
        }
    }
}

void C_KronSums_center_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *centerx, 
                               double *centery, double *PQ_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = (y[qN + i] - centery[q]) * weights[i];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += (x[p * N + i] - centerx[p]) * tmp;
        }
    }
}

void C_KronSums_center_subset(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *PQ_ans) 
{
    int qP, qN, pN;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            for (int i = 0; i < Nsubset; i++)
                 PQ_ans[qP + p] += (y[qN + subsety[i]] - centery[q]) * 
                                   (x[pN + subsetx[i]] - centerx[p]);
        }
    }
}

void C_KronSums_center_weights_subset(double *x, int N, int P, double *y, 
                                      int Q, int *weights, int *subset, 
                                      int Nsubset, double *centerx, 
                                      double *centery, double *PQ_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) PQ_ans[qP + p] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = (y[qN + subset[i]] - centery[q]) * weights[subset[i]];
             for (int p = 0; p < P; p++)
                 PQ_ans[qP + p] += (x[p * N + subset[i]] - centerx[p]) * tmp;
        }
    }
}
