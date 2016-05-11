
#include "libcoin.h"

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
    for (int i = 0; i < N; i++) ans += weights[i];
    return(ans);
}

/* sum(weights[subset]) */
int C_sum_subset(int *weights, int N, int *subset, int Nsubset)
{
    int ans = 0;
    for (int i = 0; i < Nsubset; i++) ans += weights[subset[i]];
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

/* rowSums(x) */
void C_rowSums(double *x, int N, int P, double *N_ans) 
{
    for (int i = 0; i < N; i++) {
        N_ans[i] = 0.0;
        for (int p = 0; p < P; p++)
            N_ans[i] += x[p * N + i];
    }
}

/* colSums(x) */
void C_colSums_i(int *x, int N, int P, int *P_ans) 
{
    int pN;
    
    for (int p = 0; p < P; p++) {
        P_ans[p] = 0.;
        pN = p * N;
        for (int i = 0; i < N; i++)
            P_ans[p] += x[pN + i];
    }
}

/* rowSums(x) */
void C_rowSums_i(int *x, int N, int P, int *N_ans) 
{
    for (int i = 0; i < N; i++) {
        N_ans[i] = 0;
        for (int p = 0; p < P; p++)
            N_ans[i] += x[p * N + i];
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
    int qP, qN, pN, qPp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            qPp = qP + p;
            for (int i = 0; i < Nsubset; i++)
                PQ_ans[qPp] += y[qN + subsety[i]] * x[pN + subsetx[i]];
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

/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym (double *x, int N, int P, double *PP_sym_ans) 
{
    int pN, qP, qN, SpqP;
    
    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] +=  x[qN + i] * x[pN + i];
        }
    }
}

/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_sym_weights(double *x, int N, int P, int *weights, double *PP_sym_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = x[qN + i] * weights[i];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += x[p * N + i] * tmp;
        }
    }
}

void C_KronSums_sym_subset(double *x, int N, int P, 
                           int *subset, int Nsubset, double *PP_sym_ans) 
{
    int qP, qN, pN, SpqP;
        
    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < Nsubset; i++)
                PP_sym_ans[SpqP] += x[qN + subset[i]] * x[pN + subset[i]];
        }
    }
}

void C_KronSums_sym_weights_subset(double *x, int N, int P, 
                               int *weights, int *subset, int Nsubset, 
                               double *PP_sym_ans) 
{
    int qP, qN;
    double tmp;
        
    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = x[qN + subset[i]] * weights[subset[i]];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += x[p * N + subset[i]] * tmp;
        }
    }
}

void C_KronSums_sym_2dweights(double *x, int N, int P, 
                          int *weights, double *PP_sym_ans) 
{
    int qPp, qN, pNi;
        
    for (int p = 0; p < P; p++) {
        for (int q = 0; q <= p; q++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            qPp = S(p, q, P);
            qN = q * N;
            for (int i = 0; i < N; i++) {
                pNi = p * N + i;
                for (int m = 0; m < N; m++)
                      PP_sym_ans[qPp] += x[qN + m] * x[pNi] * weights[m * N + i];
            }
        }
    }
}



/* sum_i (t(x[i,] - center) %*% (x[i,] - center)) */
void C_KronSums_sym_center(double *x, int N, int P, 
                       double *center, double *PP_sym_ans) 
{
    int qP, qN, pN, SpqP;
        
    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < N; i++)
                 PP_sym_ans[SpqP] += (x[qN + i] - center[q]) * (x[pN + i] - center[p]);
        }
    }
}

void C_KronSums_sym_center_weights(double *x, int N, int P, 
                               int *weights, double *center, 
                               double *PP_sym_ans) 
{
    int qP, qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < N; i++) {
             tmp = (x[qN + i] - center[q]) * weights[i];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += (x[p * N + i] - center[p]) * tmp;
        }
    }
}

void C_KronSums_sym_center_subset(double *x, int N, int P, 
                              int *subset, int Nsubset, 
                              double *center, double *PP_sym_ans) 
{
    int qP, qN, pN, SpqP;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) {
            PP_sym_ans[S(p, q, P)] = 0.0;
            pN = p * N;
            SpqP = S(p, q, P);
            for (int i = 0; i < Nsubset; i++)
                PP_sym_ans[SpqP] += (x[qN + subset[i]] - center[q]) *
                                    (x[pN + subset[i]] - center[p]); 
                                                   
        }
    }
}

void C_KronSums_sym_center_weights_subset(double *x, int N, int P, 
                                      int *weights, int *subset, 
                                      int Nsubset, double *center,
                                      double *PP_sym_ans) 
{
    int qP, qN;
    double tmp;

    for (int q = 0; q < P; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p <= q; p++) PP_sym_ans[S(p, q, P)] = 0.0;
        for (int i = 0; i < Nsubset; i++) {
             tmp = (x[qN + subset[i]] - center[q]) * weights[subset[i]];
             for (int p = 0; p <= q; p++)
                 PP_sym_ans[S(p, q, P)] += (x[p * N + subset[i]] - center[p]) * tmp;
        }
    }
}

/* sum_i (t(x[i,] - centerx) %*% (y[i,] - centery)) */
void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *PQ_ans) 
{
    int qP, qN, pN, qPp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            qPp = qP + p;
            for (int i = 0; i < N; i++)
                 PQ_ans[qPp] += (y[qN + i] - centery[q]) * 
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
    int qP, qN, pN, qPp;
        
    for (int q = 0; q < Q; q++) {
        qN = q * N;
        qP = q * P;
        for (int p = 0; p < P; p++) {
            PQ_ans[qP + p] = 0.0;
            pN = p * N;
            qPp = qP + p;
            for (int i = 0; i < Nsubset; i++)
                 PQ_ans[qPp] += (y[qN + subsety[i]] - centery[q]) * 
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

/* tapply(1:nrow(y), ix, function(i) colSums(y[i,])) */
void C_tapplySum(double *y, int N, int Q, int *ix, int Nx, double *NxQ_ans)
{
    int iN, ixi;
   
    for (int q = 0; q < Nx * Q; q++) NxQ_ans[q] = 0.0;
   
    for (int i = 0; i < N; i++) {
        iN = i * N;
        ixi = ix[i] - 1; /* ix[i] == 0 means NA */
        if (ixi >= 0) {
            for (int q = 0; q < Q; q++)
                NxQ_ans[q * Nx + ixi] += y[iN + q];
        }
    }
}

void C_tapplySum_weights(double *y, int N, int Q, int *ix, int Nx, 
                         int *weights, double *NxQ_ans)
{
    int iN, ixi, wi;
   
    for (int q = 0; q < Nx * Q; q++) NxQ_ans[q] = 0.0;
   
    for (int i = 0; i < N; i++) {
        iN = i * N;
        wi = weights[i];
        ixi = ix[i] - 1;
        if (ixi >= 0) {
            for (int q = 0; q < Q; q++)
                NxQ_ans[q * Nx + ixi] += wi * y[iN + q];
        }
    }
}

void C_tapplySum_subset(double *y, int N, int Q, int *ix, int Nx, 
                        int *subset, int Nsubset, double *NxQ_ans)
{
    int iN, ixi;
   
    for (int q = 0; q < Nx * Q; q++) NxQ_ans[q] = 0.0;
   
    for (int i = 0; i < Nsubset; i++) {
        iN = subset[i] * N;
        ixi = ix[subset[i]] - 1;
        if (ixi >= 0) {
            for (int q = 0; q < Q; q++)
                NxQ_ans[q * Nx + ixi] += y[iN + q];
        }
    }
}

void C_tapplySum_weights_subset(double *y, int N, int Q, int *ix, int Nx, 
                                int *weights, int *subset, int Nsubset, 
                                double *NxQ_ans)
{
    int iN, wi, ixi;
   
    for (int q = 0; q < Nx * Q; q++) NxQ_ans[q] = 0.0;
   
    for (int i = 0; i < Nsubset; i++) {
        iN = subset[i] * N;
        wi = weights[subset[i]];
        ixi = ix[subset[i]] - 1;
        if (ixi >= 0) {
            for (int q = 0; q < Q; q++)
                NxQ_ans[q * Nx + ixi] += wi * y[iN + q];
        }
    }
}

void C_tapplySum_2d(double *y, int M, int Q, int N, 
                    int *weights, double *Nm1Q_ans)
{
    for (int q = 0; q < (N - 1) * Q; q++) Nm1Q_ans[q] = 0.0;
   
    for (int m = 1; m < M; m++) { /* m = 0 means NA */
        for (int i = 1; i < N; i++) { /* i = 0 means NA */
            for (int q = 0; q < Q; q++) 
                Nm1Q_ans[q * (N - 1) + (i - 1)] += 
                    weights[q * M + i] * y[q * M + m];
        }
    }
}
