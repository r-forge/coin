
#include "libcoin_internal.h"
#include "Utils.h"

/* Variables:
  ix:		integer vector of length N with elements 0...(Lx - 1)
  iy:		integer vector of length N with elements 0...(Ly - 1)
  weights:	integer vector of length N
  subset:       an integer Nsubset vector with elements 0...(N - 1)
  block:	an integer N vector with elements 1...Lb
  LxLy_ans:	return value, integer matrix Lx x Ly 
  LxLyLb_ans:	return value, integer array Lx x Ly x Lb 
  Lx_ans:	return value, integer vector of length Lx
*/


/* table(ix) */
void C_1dtable_
(
    int *ix, 
    int Lx, 
    int N, 
    int *Lx_ans
) {

    for (int i = 0; i < Lx; i++) Lx_ans[i] = 0;

    for (int i = 0; i < N; i++) Lx_ans[ix[i]]++;
}

/* table(ix[subset]) */
void C_1dtable_subset
(
    int *ix, 
    int Lx, 
    int *subset,
    int Nsubset, 
    int *Lx_ans
) {

    for (int i = 0; i < Lx; i++) Lx_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
        Lx_ans[ix[subset[i]]]++;
}

/* xtabs(weights ~ ix) */
void C_1dtable_weights
(
    int *ix, 
    int Lx, 
    int *weights, 
    int N, 
    int *Lx_ans
) {

    for (int i = 0; i < Lx; i++) Lx_ans[i] = 0;
    
    for (int i = 0; i < N; i++)  
        Lx_ans[ix[i]] += weights[i];
}

/* xtabs(weights ~ ix, subset = subset) */
void C_1dtable_weights_subset
(
    int *ix, 
    int Lx, 
    int *weights, 
    int *subset, 
    int Nsubset,
    int *Lx_ans
) {

    for (int i = 0; i < Lx; i++) Lx_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
          Lx_ans[ix[subset[i]]] += weights[subset[i]];
}

/* table(ix, iy) */
void C_2dtable_
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int N, 
    int *LxLy_ans
) {

    for (int i = 0; i < Lx * Ly; i++) LxLy_ans[i] = 0;
    
    for (int i = 0; i < N; i++)  
        LxLy_ans[ix[i] + iy[i] * Lx]++;

}

/* table(ix[subset], iy[subset]) */
void C_2dtable_subset
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *subset, 
    int Nsubset, 
    int *LxLy_ans
) {

    for (int i = 0; i < Lx * Ly; i++) LxLy_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
        LxLy_ans[ix[subset[i]] + iy[subset[i]] * Lx]++;

}

/* xtabs(weights ~ ix + iy) */
void C_2dtable_weights
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *weights,
    int N, 
    int *LxLy_ans
) {

    for (int i = 0; i < Lx * Ly; i++) LxLy_ans[i] = 0;
    
    for (int i = 0; i < N; i++)
        LxLy_ans[ix[i] + iy[i] * Lx] += weights[i];
}

/* xtabs(weights ~ ix + iy, subset = subset) */
void C_2dtable_weights_subset
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *weights,
    int *subset, 
    int Nsubset, 
    int *LxLy_ans
) {

    for (int i = 0; i < Lx * Ly; i++) LxLy_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)
         LxLy_ans[ix[subset[i]] + iy[subset[i]] * Lx] += weights[subset[i]];
}

/* table(ix, iy, block) w/o NAs in block, ie block > 0 */
void C_2dtable_block
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *block, 
    int Lb, 
    int N, 
    int *LxLyLb_ans
) {

    int LxLy = Lx * Ly;

    for (int i = 0; i < LxLy * Lb; i++) LxLyLb_ans[i] = 0;
    
    for (int i = 0; i < N; i++)  
        LxLyLb_ans[(block[i] - 1) * LxLy + ix[i] + iy[i] * Lx]++;

}

/* table(ix[subset], iy[subset], block[subset]) w/o NAs in block, ie block > 0 */
void C_2dtable_subset_block
(
    int *ix, 
    int Lx,
    int *iy, 
    int Ly, 
    int *subset, 
    int Nsubset, 
    int *block, 
    int Lb, 
    int *LxLyLb_ans
) {

    int LxLy = Lx * Ly;

    for (int i = 0; i < LxLy * Lb; i++) LxLyLb_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)  
        LxLyLb_ans[(block[subset[i]] - 1) * LxLy + 
                   ix[subset[i]] + iy[subset[i]] * Lx]++;

}

/* xtabs(weights ~ ix + iy + block) w/o NAs in block, ie block > 0 */
void C_2dtable_weights_block
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *weights,
    int *block, 
    int Lb, 
    int N, 
    int *LxLyLb_ans
) {

    int LxLy = Lx * Ly;

    for (int i = 0; i < LxLy * Lb; i++) LxLyLb_ans[i] = 0;
    
    for (int i = 0; i < N; i++)
        LxLyLb_ans[(block[i] - 1) * LxLy + ix[i] + iy[i] * Lx] += weights[i];
}

/* xtabs(weights ~ ix + iy + block, subset = subset) w/o NAs in block, ie block > 0 */
void C_2dtable_weights_subset_block
(
    int *ix, 
    int Lx, 
    int *iy, 
    int Ly, 
    int *weights, 
    int *subset, 
    int Nsubset, 
    int *block, 
    int Lb, 
    int *LxLyLb_ans
) {

    int LxLy = Lx * Ly;

    for (int i = 0; i < LxLy * Lb; i++) LxLyLb_ans[i] = 0;
    
    for (int i = 0; i < Nsubset; i++)
         LxLyLb_ans[(block[subset[i]] - 1) * LxLy + 
                    ix[subset[i]] + iy[subset[i]] * Lx] += weights[subset[i]];
}

void RC_2dtable
(
    SEXP ix, 
    SEXP iy, 
    SEXP weights, 
    SEXP subset, 
    SEXP block, 
    int *LxLyLb_ans
) {

    if (LENGTH(block) == 0) {
        if ((LENGTH(weights) == 0) && (LENGTH(subset) == 0))
            C_2dtable_(INTEGER(ix), NLEVELS(ix) + 1, 
                       INTEGER(iy), NLEVELS(iy) + 1, 
                       LENGTH(ix), LxLyLb_ans);
        if ((LENGTH(weights) > 0) && (LENGTH(subset) == 0))
            C_2dtable_weights(INTEGER(ix), NLEVELS(ix) + 1, 
                              INTEGER(iy), NLEVELS(iy) + 1, 
                              INTEGER(weights), LENGTH(weights), LxLyLb_ans);
        if ((LENGTH(weights) == 0) && (LENGTH(subset) > 0))
            C_2dtable_subset(INTEGER(ix), NLEVELS(ix) + 1, 
                             INTEGER(iy), NLEVELS(iy) + 1, 
                             INTEGER(subset), LENGTH(subset), LxLyLb_ans);
        if ((LENGTH(weights) > 0) && (LENGTH(subset) > 0))
            C_2dtable_weights_subset(INTEGER(ix), NLEVELS(ix) + 1, 
                                     INTEGER(iy), NLEVELS(iy) + 1, 
                                     INTEGER(weights), INTEGER(subset), 
                                     LENGTH(subset), LxLyLb_ans);
    } else {
        if ((LENGTH(weights) == 0) && (LENGTH(subset) == 0))
            C_2dtable_block(INTEGER(ix), NLEVELS(ix) + 1, 
                            INTEGER(iy), NLEVELS(iy) + 1, 
                            INTEGER(block), NLEVELS(block), LENGTH(block), 
                            LxLyLb_ans);
        if ((LENGTH(weights) > 0) && (LENGTH(subset) == 0))
            C_2dtable_weights_block(INTEGER(ix), NLEVELS(ix) + 1, 
                                    INTEGER(iy), NLEVELS(iy) + 1, 
                                    INTEGER(weights),  
                                    INTEGER(block), NLEVELS(block), 
                                    LENGTH(weights), LxLyLb_ans);
        if ((LENGTH(weights) == 0) && (LENGTH(subset) > 0))
            C_2dtable_subset_block(INTEGER(ix), NLEVELS(ix) + 1, 
                                   INTEGER(iy), NLEVELS(iy) + 1, 
                                   INTEGER(subset), LENGTH(subset), 
                                   INTEGER(block), NLEVELS(block), 
                                   LxLyLb_ans);
        if ((LENGTH(weights) > 0) && (LENGTH(subset) > 0))
            C_2dtable_weights_subset_block(INTEGER(ix), NLEVELS(ix) + 1, 
                                           INTEGER(iy), NLEVELS(iy) + 1, 
                                           INTEGER(weights), INTEGER(subset),
                                           LENGTH(subset), INTEGER(block), 
                                           NLEVELS(block), LxLyLb_ans);
    }
}
