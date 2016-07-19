
/* Variables
   x:           a double or integer N x P or Lx x P matrix
                the first row is == 0 in the latter case
   y:           a double N x Q or Ly x Q matrix
                the first two is == 0 in the latter case
   ix:          an integer vector of length N with elements 0...(Lx - 1)
                0 means NA
   weights:     an integer N vector 
   weights2d:   an integer Lx x Ly vector 
   subset:      an integer Nsubset vector with elements 0...(N - 1)
   subsety:     an integer Nsubset vector with elements 0...(N - 1)
   subsety:     an integer Nsubset vector with elements 0...(N - 1)
   centerx:	a double P vector centering the columns of x
   centery:	a double Q vector centering the columns of y
   PQ_ans:      return value, a double P x Q matrix
   P_ans:       return value, a double or integer P vector
   N_ans:       return value, a double or integer N vector
   Q_ans:       return value, a double Q vector
   PP_sym_ans:	return value, a symmetric double P x P matrix in lower packed format
   PQ_ans:	return value, a double P x Q matrix 
   LxQ_ans:	return value, a double Lx x Q matrix
   Lx1Q_ans:	return value, a double (Lx - 1) x Q matrix
*/
                                                   

/* sum(weights) */
int C_sum_(int *weights, int N);
/* sum(weights[subset]) */
int C_sum_subset(int *weights, int N, int *subset, int Nsubset);
/* colSums(x) */
void C_colSums_(double *x, int N, int P, double *P_ans);
/* rowSums(x) */
void C_rowSums_(double *x, int N, int P, double *N_ans); 
/* colSums(x) integer version */
void C_colSums_i(int *x, int N, int P, int *P_ans);
/* rowSums(x) integer version */
void C_rowSums_i(int *x, int N, int P, int *N_ans); 
/* colSums(x * weights) */
void C_colSums_weights(double *x, int N, int P, int *weights, double *P_ans); 
/* colSums(x[subsetx,]) */
void C_colSums_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                      double *P_ans); 
/* colSums(x[subsetx,] * weights[subsetx]) */
void C_colSums_weights_subset(double *x, int N, int P, int *weights, 
                              int *subsetx, int Nsubset, double *P_ans); 
/* colSums(x^2) */
void C_colSums2_(double *x, int N, int P, double *P_ans); 
/* colSums(x^2 * weights) */
void C_colSums2_weights(double *x, int N, int P, int *weights, double *P_ans); 
/* colSums(x[subsetx,]^2) */
void C_colSums2_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                       double *P_ans); 
/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_weights_subset(double *x, int N, int P, int *weights, 
                               int *subsetx, int Nsubset, double *P_ans); 
/* colSums((x-center)^2) */
void C_colSums2_center_(double *x, int N, int P, double *centerx, double *P_ans); 
/* colSums((x-center)^2 * weights) */
void C_colSums2_center_weights(double *x, int N, int P, int *weights, 
                               double *centerx, double *P_ans); 
/* colSums((x[subsetx,] - center)^2) */
void C_colSums2_center_subset(double *x, int N, int P, int *subsetx, 
                              int Nsubset, double *centerx, double *P_ans); 
/* colSums((x[subsetx,]-center)^2 * weights[subsetx]) */
void C_colSums2_center_weights_subset(double *x, int N, int P, int *weights, 
                                      int *subsetx, int Nsubset, 
                                      double *centerx, double *P_ans); 
/* sum_i (t(x[i,]) %*% y[i,]) */
void C_KronSums_(double *x, int N, int P, double *y, int Q, double *PQ_ans); 
/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        int *weights, double *PQ_ans); 
/* sum_i (t(x[subsetx[i],]) %*% y[subsety[i],]) */
void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, 
                       double *PQ_ans); 
/* sum_i weights[subset[i]] (t(x[subset[i],]) %*% y[subset[i],]) */
void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               int *weights, int *subset, int Nsubset, 
                               double *PQ_ans); 
/* sum_i,j weights2d[i, j] * t(x[i,]) %*% y[j,]) */
void C_KronSums_2dweights(double *x, int Lx, int P, double *y, int Ly, int Q, 
                          int *weights2d, double *PQ_ans); 
/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_(double *x, int N, int P, double *PP_sym_ans);
/* sum_i weights[i] * (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym_weights(double *x, int N, int P, int *weights, 
                            double *PP_sym_ans); 
/* sum_i (t(x[subset[i],]) %*% x[subset[i],]) */
void C_KronSums_sym_subset(double *x, int N, int P, 
                           int *subset, int Nsubset, double *PP_sym_ans); 
/* sum_i weights[subset[i]] (t(x[subsetx[i],]) %*% y[subsety[i],]) */
void C_KronSums_sym_weights_subset(double *x, int N, int P, 
                               int *weights, int *subset, int Nsubset, 
                               double *PP_sym_ans); 
/* sum_i (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
void C_KronSums_sym_center_(double *x, int N, int P, 
                            double *centerx, double *PP_sym_ans); 
/* sum_i weights[i] (t(x[i,] - centerx) %*% (x[i,] - centerx)) */
void C_KronSums_sym_center_weights(double *x, int N, int P, 
                                   int *weights, double *centerx, 
                                   double *PP_sym_ans); 
/* sum_i (t(x[subset[i],] - centerx) %*% (x[subset[i],] - centerx)) */
void C_KronSums_sym_center_subset(double *x, int N, int P, 
                                  int *subset, int Nsubset, 
                                  double *centerx, double *PP_sym_ans); 
/* sum_i weights[subset[i]] (t(x[subset[i],] - centerx) %*% 
                            (x[subset[i],] - centerx)) */
void C_KronSums_sym_center_weights_subset(double *x, int N, int P, 
                                          int *weights, int *subset, 
                                          int Nsubset, double *centerx,
                                          double *PP_sym_ans); 
/* tapply(1:nrow(y), ix, function(i) colSums(y[i,])) */
void C_tapplySum_(double *y, int N, int Q, int *ix, int Lx, double *LxQ_ans);
/* tapply(1:nrow(y), ix, function(i) colSums(weights[i] * y[i,])) */
void C_tapplySum_weights(double *y, int N, int Q, int *ix, int Lx, 
                         int *weights, double *LxQ_ans);
/* tapply((1:nrow(y))[subsety], ix[subsetx], 
          function(i) colSums(y[i,])) */
void C_tapplySum_subset(double *y, int N, int Q, int *ix, int Lx, 
                        int *subsetx, int *subsety, int Nsubset, double *LxQ_ans);
/* tapply((1:nrow(y))[subset], ix[subset], 
          function(i) colSums(weights[i] * y[i,])) */
void C_tapplySum_weights_subset(double *y, int N, int Q, int *ix, int Lx, 
                                int *weights, int *subset, int Nsubset, 
                                double *LxQ_ans);
/* sum(weights[i, ] * y[, q]) forall i = 1, ... Lx and q = 0, ..., Q */
void C_tapplySum_2d(double *y, int Ly, int Q, int N, 
                    int *weights2d, double *Lx1Q_ans);
