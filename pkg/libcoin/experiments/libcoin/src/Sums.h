
void C_sum(double *x, int N, double *ans);
void C_sum_subset(double *x, int N, int *subset, int Nsubset, double *ans);
void C_colSums(double *x, int N, int P, double *ans);
void C_colSums_weights(double *x, int N, int P, double *weights, double *ans);
void C_colSums_subset(double *x, int N, int P, int *subset, int Nsubset, double *ans);
void C_colSums_weights_subset(double *x, int N, int P, double *weights, 
                              int *subset, int Nsubset, double *ans);
void C_KronSums(double *x, int N, int P, double *y, int Q, double *ans);
void C_KronSums_weights(double *x, int N, int P, double *y, int Q, 
                        double *weights, double *ans);
void C_KronSums_subset(double *x, int N, int P, double *y, int Q, 
                       int *subsetx, int *subsety, int Nsubset, double *ans);
void C_KronSums_weights_subset(double *x, int N, int P, double *y, int Q, 
                               double *weights, int *subset, int Nsubset, 
                               double *ans);
void C_KronSums_2dweights(double *x, int N, int P, double *y, int M, int Q, 
                          double *weights, double *ans);
void C_KronSums_center(double *x, int N, int P, double *y, int Q, 
                       double *centerx, double *centery, double *ans);
void C_KronSums_weights_center(double *x, int N, int P, double *y, int Q, 
                               double *weights, double *centerx, 
                               double *centery, double *ans);
void C_KronSums_subset_center(double *x, int N, int P, double *y, int Q, 
                              int *subsetx, int *subsety, int Nsubset, 
                              double *centerx, double *centery, double *ans);
void C_KronSums_weights_subset_center(double *x, int N, int P, double *y, 
                                      int Q, double *weights, 
                                      int *subset, int Nsubset, 
                                      double *centerx, double *centery, double *ans);
