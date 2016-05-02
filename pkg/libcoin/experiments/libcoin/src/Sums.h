
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
int C_sum(int *weights, int N);
/* sum(weights[subset]) */
int C_sum_subset(int *weights, int N, int *subset, int Nsubset);
/* colSums(x) */
void C_colSums(double *x, int N, int P, double *P_ans);
/* colSums(x * w) */
void C_colSums_weights(double *x, int N, int P, int *weights, double *P_ans);
/* colSums(x[subsetx,]) */
void C_colSums_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                      double *P_ans);
/* colSums(x[subsetx,] * weights[subsetx]) */
void C_colSums_weights_subset(double *x, int N, int P, int *weights, 
                              int *subsetx, int Nsubset, double *P_ans); 
/* colSums(x^2) */
void C_colSums2(double *x, int N, int P, double *P_ans);
/* colSums(x^2 * w) */
void C_colSums2_weights(double *x, int N, int P, int *weights, double *P_ans); 
/* colSums(x[subsetx,]^2) */
void C_colSums2_subset(double *x, int N, int P, int *subsetx, int Nsubset, 
                       double *P_ans); 
/* colSums(x[subsetx,]^2 * weights[subsetx]) */
void C_colSums2_weights_subset(double *x, int N, int P, int *weights, 
                               int *subsetx, int Nsubset, double *P_ans); 
/* colSums((x-center)^2) */
void C_colSums2_center(double *x, int N, int P, double *centerx, double *P_ans); 
/* colSums((x-center)^2 * w) */
void C_colSums2_center_weights(double *x, int N, int P, int *weights, 
                               double *centerx, double *P_ans); 
/* colSums((x[subsetx,] - center)^2) */
void C_colSums2_center_subset(double *x, int N, int P, int *subsetx, 
                              int Nsubset, double *centerx, double *P_ans); 
/* colSums((x[subsetx,]-center)^2 * weights[subsetx]) */
void C_colSums2_center_weights_subset(double *x, int N, int P, int *weights, 
                                      int *subsetx, int Nsubset, double *centerx, 
                                      double *P_ans); 
/* sum_i (t(x[i,]) %*% y[i,]) */
void C_KronSums(double *x, int N, int P, double *y, int Q, double *PQ_ans); 
/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        int *weights, double *PQ_ans); 
void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, double *PQ_ans); 
void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               int *weights, int *subset, int Nsubset, 
                               double *PQ_ans); 
void C_KronSums_2dweights(double *x, int N, int P, double *y, int M, int Q, 
                          int *weights, double *PQ_ans); 
/* sum_i (t(x[i,]) %*% x[i,]) */
void C_KronSums_sym(double *x, int N, int P, double *PP_sym_ans); 
/* sum_i weights[i] * (t(x[i,]) %*% y[i,]) */
void C_KronSums_sym_weights(double *x, int N, int P, 
                        int *weights, double *PP_sym_ans); 
void C_KronSums_sym_subset(double *x, int N, int P, 
                       int *subset, int Nsubset, double *PP_sym_ans); 
void C_KronSums_sym_weights_subset(double *x, int N, int P, 
                               int *weights, int *subset, int Nsubset, 
                               double *PP_sym_ans); 
void C_KronSums_sym_2dweights(double *x, int N, int P, 
                          int *weights, double *PP_sym_ans); 
/* sum_i (t(x[i,] - center) %*% (x[i,] - center)) */
void C_KronSums_sym_center(double *x, int N, int P, 
                       double *center, double *PP_sym_ans); 
void C_KronSums_sym_center_weights(double *x, int N, int P, 
                               int *weights, double *center, 
                               double *PP_sym_ans); 
void C_KronSums_sym_center_subset(double *x, int N, int P, 
                              int *subsetx, int Nsubset, 
                              double *center, double *PP_sym_ans); 
void C_KronSums_sym_center_weights_subset(double *x, int N, int P, 
                                      int *weights, int *subset, 
                                      int Nsubset, double *center,
                                      double *PP_sym_ans); 
/* sum_i (t(x[i,] - centerx) %*% (y[i,] - centery)) */
void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *PQ_ans); 
void C_KronSums_center_weights(double *x, int N, int P, double *y, int Q, 
                               int *weights, double *centerx, 
                               double *centery, double *PQ_ans); 
void C_KronSums_center_subset(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *PQ_ans); 
void C_KronSums_center_weights_subset(double *x, int N, int P, double *y, 
                                      int Q, int *weights, int *subset, 
                                      int Nsubset, double *centerx, 
                                      double *centery, double *PQ_ans); 
/* tapply(1:nrow(y), ix, function(i) colSums(y[i,])) */
void C_tapplySum(double *y, int N, int Q, int *ix, int Nx, double *NxQ_ans);
void C_tapplySum_weights(double *y, int N, int Q, int *ix, int Nx, 
                         int *weights, double *NxQ_ans);
void C_tapplySum_subset(double *y, int N, int Q, int *ix, int Nx, 
                        int *subset, int Nsubset, double *NxQ_ans);
void C_tapplySum_weights_subset(double *y, int N, int Q, int *ix, int Nx, 
                                int *weights, int *subset, int Nsubset, 
                                double *NxQ_ans);
void C_tapplySum_2d(double *y, int M, int Q, int Nx, 
                    int *weights, double *NxQ_ans);
