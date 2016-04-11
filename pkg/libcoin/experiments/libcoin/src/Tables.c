
#include "libcoin.h"
#include "helpers.h"

void C_2dtable_set_zero(int Nx, int Ny, int *NxNy_ans) {

    for (int j = 0; j < Nx; j++) {
         for (int k = 0; k < Ny; k++)
              NxNy_ans[k * Nx + j] = 0;
    }
}

/* table(ix, iy) */
void C_2dtable(int *ix, int Nx, int *iy, int Ny, int N, int *NxNy_ans) {

    C_2dtable_set_zero(Nx, Ny, NxNy_ans);
    
    for (int i = 0; i < N; i++)  
        NxNy_ans[ix[i] + iy[i] * Nx]++;

}

/* table(ix[subset], iy[subset]) */
void C_2dtable_subset(int *ix, int Nx, int *iy, int Ny, int *subset, 
                      int Nsubset, int *NxNy_ans) {

    C_2dtable_set_zero(Nx, Ny, NxNy_ans);
    
    for (int i = 0; i < Nsubset; i++)  
        NxNy_ans[ix[subset[i]] + iy[subset[i]] * Nx]++;

}

/* xtabs(weights ~ ix + iy) */
void C_2dtable_weights(int *ix, int Nx, int *iy, int Ny, int *weights,
                       int N, int *NxNy_ans) {

    C_2dtable_set_zero(Nx, Ny, NxNy_ans);
    
    for (int i = 0; i < N; i++)
        NxNy_ans[ix[i] + iy[i] * Nx] += weights[i];
}

/* xtabs(weights ~ ix + iy, subset = subset) */

void C_2dtable_weights_subset(int *ix, int Nx, int *iy, int Ny, int *weights,
                              int *subset, int Nsubset, int *NxNy_ans) {

    C_2dtable_set_zero(Nx, Ny, NxNy_ans);
    
    for (int i = 0; i < Nsubset; i++)
         NxNy_ans[ix[subset[i]] + iy[subset[i]] * Nx] += weights[subset[i]];
}

void C_1dtable_set_zero(int Ny, int *Ny_ans) {
    for (int i = 0; i < Ny; i++) Ny_ans[i] = 0;
}

void C_1dtable(int *iy, int Ny, int N, int *Ny_ans) {

    C_1dtable_set_zero(Ny, Ny_ans);
    for (int i = 0; i < N; i++) Ny_ans[iy[i]]++;
}

void C_1dtable_subset(int *iy, int Ny, int *subset, int Nsubset, int *Ny_ans) {

    C_1dtable_set_zero(Ny, Ny_ans);
    
    for (int i = 0; i < Nsubset; i++)  
        Ny_ans[iy[subset[i]]]++;
}

void C_1dtable_weights(int *iy, int Ny, int *weights, int N, int *Ny_ans) {

    C_1dtable_set_zero(Ny, Ny_ans);
    
    for (int i = 0; i < N; i++)  
        Ny_ans[iy[i]] += weights[i];
}

void C_1dtable_weights_subset(int *iy, int Ny, int *weights, int *subset, 
                              int Nsubset, int *Ny_ans) {

    C_1dtable_set_zero(Ny, Ny_ans);
    
    for (int i = 0; i < Nsubset; i++)  
          Ny_ans[iy[subset[i]]] += weights[subset[i]];
          
}
