
#include "libcoin.h"
#include "helpers.h"

/* Variables:
  ix:		integer vector of length N with elements 0...Nx
  iy:		integer vector of length N with elements 0...Ny
  weights:	integer vector of length N
  subset:	integer vector of length N
  NxNy_ans:	integer matrix Nx x Ny 
  Nx_ans:	integer vector of length Nx
*/

/* table(ix, iy) */
void C_2dtable(int *ix, int Nx, int *iy, int Ny, int N, int *NxNy_ans)
{
    for (int i = 0; i < Nx * Ny; i++) NxNy_ans[i] = 0;
    
    for (int i = 0; i < N; i++)  
        NxNy_ans[ix[i] + iy[i] * Nx]++;

}

/* table(ix[subset], iy[subset]) */
void C_2dtable_subset(int *ix, int Nx, int *iy, int Ny, int *subset, 
                      int Nsubset, int *NxNy_ans) 
{
    for (int i = 0; i < Nx * Ny; i++) NxNy_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
        NxNy_ans[ix[subset[i]] + iy[subset[i]] * Nx]++;

}

/* xtabs(weights ~ ix + iy) */
void C_2dtable_weights(int *ix, int Nx, int *iy, int Ny, int *weights,
                       int N, int *NxNy_ans) 
{
    for (int i = 0; i < Nx * Ny; i++) NxNy_ans[i] = 0;
    
    for (int i = 0; i < N; i++)
        NxNy_ans[ix[i] + iy[i] * Nx] += weights[i];
}

/* xtabs(weights ~ ix + iy, subset = subset) */
void C_2dtable_weights_subset(int *ix, int Nx, int *iy, int Ny, int *weights,
                              int *subset, int Nsubset, int *NxNy_ans) 
{
    for (int i = 0; i < Nx * Ny; i++) NxNy_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)
         NxNy_ans[ix[subset[i]] + iy[subset[i]] * Nx] += weights[subset[i]];
}

/* table(ix) */
void C_1dtable(int *ix, int Nx, int N, int *Nx_ans) 
{
    for (int i = 0; i < Nx; i++) Nx_ans[i] = 0;

    for (int i = 0; i < N; i++) Nx_ans[ix[i]]++;
}

/* table(ix[subset]) */
void C_1dtable_subset(int *ix, int Nx, int *subset, int Nsubset, int *Nx_ans) 
{
    for (int i = 0; i < Nx; i++) Nx_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
        Nx_ans[ix[subset[i]]]++;
}

/* xtabs(weights ~ ix) */
void C_1dtable_weights(int *ix, int Nx, int *weights, int N, int *Nx_ans) 
{
    for (int i = 0; i < Nx; i++) Nx_ans[i] = 0;
    
    for (int i = 0; i < N; i++)  
        Nx_ans[ix[i]] += weights[i];
}

/* xtabs(weights ~ ix, subset = subset) */
void C_1dtable_weights_subset(int *ix, int Nx, int *weights, int *subset, 
                              int Nsubset, int *Nx_ans) 
{
    for (int i = 0; i < Nx; i++) Nx_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
          Nx_ans[ix[subset[i]]] += weights[subset[i]];
}
