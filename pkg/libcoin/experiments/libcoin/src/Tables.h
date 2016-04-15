
/* Variables:
  ix:		integer vector of length N with elements 0...Nx
  iy:		integer vector of length N with elements 0...Ny
  weights:	integer vector of length N
  subset:	integer vector of length N
  NxNy_ans:	integer matrix Nx x Ny 
  Nx_ans:	integer vector of length Nx
*/

/* table(ix, iy) */
void C_2dtable(int *ix, int Nx, int *iy, int Ny, int N, int *NxNy_ans);
/* table(ix[subset], iy[subset]) */
void C_2dtable_subset(int *ix, int Nx, int *iy, int Ny, int *subset, 
                      int Nsubset, int *NxNy_ans);
/* xtabs(weights ~ ix + iy) */
void C_2dtable_weights(int *ix, int Nx, int *iy, int Ny, int *weights,
                       int N, int *NxNy_ans);
/* xtabs(weights ~ ix + iy, subset = subset) */
void C_2dtable_weights_subset(int *ix, int Nx, int *iy, int Ny, int *weights,
                              int *subset, int Nsubset, int *NxNy_ans);
/* table(ix) */
void C_1dtable(int *ix, int Nx, int N, int *Nx_ans);
/* table(ix[subset]) */
void C_1dtable_subset(int *ix, int Nx, int *subset, int Nsubset, int *Nx_ans);
/* xtabs(weights ~ ix) */
void C_1dtable_weights(int *ix, int Nx, int *weights, int N, int *Nx_ans);
/* xtabs(weights ~ ix, subset = subset) */
void C_1dtable_weights_subset(int *ix, int Nx, int *weights, int *subset, 
                              int Nsubset, int *Nx_ans);

void C_2dtable_block(int *ix, int Nx, int *iy, int Ny, int *block, int Nlevels, int N, int *NxNyNlevels_ans);
void C_2dtable_subset_block(int *ix, int Nx, int *iy, int Ny, int *subset, int Nsubset, int *block, int Nlevels,
                            int *NxNy_ans);
void C_2dtable_weights_block(int *ix, int Nx, int *iy, int Ny, int *weights,int *block, int Nlevels,
                       int N, int *NxNy_ans);
void C_2dtable_weights_subset_block(int *ix, int Nx, int *iy, int Ny, int *weights,
                              int *subset, int Nsubset, int *block, int Nlevels, int *NxNy_ans);
