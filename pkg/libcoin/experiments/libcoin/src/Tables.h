
void C_2dtable(int *ix, int Nx, int *iy, int Ny, int N, int *NxNy_ans);
void C_2dtable_subset(int *ix, int Nx, int *iy, int Ny, int *subset, 
                      int Nsubset, int *NxNy_ans);
void C_2dtable_weights(int *ix, int Nx, int *iy, int Ny, int *weights,
                       int N, int *NxNy_ans);
void C_2dtable_weights_subset(int *ix, int Nx, int *iy, int Ny, int *weights,
                              int *subset, int Nsubset, int *NxNy_ans);
void C_1dtable(int *iy, int Ny, int N, int *Ny_ans);
void C_1dtable_subset(int *iy, int Ny, int *subset, int Nsubset, int *Ny_ans);
void C_1dtable_weights(int *iy, int Ny, int *weights, int N, int *Ny_ans);
void C_1dtable_weights_subset(int *iy, int Ny, int *weights, int *subset, 
                              int Nsubset, int *Ny_ans);
