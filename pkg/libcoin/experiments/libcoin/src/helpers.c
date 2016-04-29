
#include "libcoin.h"
#include "Tables.h"
#include "Sums.h"

int C_nlevels (SEXP x) 
{
   return(LENGTH(getAttrib(x, R_LevelsSymbol)));
}

int NROW (SEXP x) 
{
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(LENGTH(x));
    return(INTEGER(a)[0]);
}
    
int NCOL (SEXP x) 
{
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    return(INTEGER(a)[1]);
}

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans)
{
    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++)
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
            }
        }
    }
}  

void C_kronecker_sym (const double *A, const int m, 
                      const double *B, const int r, 
                      double *ans)
{
    int i, j, k, l, mr, js, ir, s, n, tmp;
    double y;

    mr = m * r;
    n = m;
    s = r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j <= i; j++) {
            js = j * s;
            y = A[S(i, j, m)]; 
            for (k = 0; k < r; k++) {
                for (l = 0; l < (j < i ? s : k + 1); l++) {
                    ans[S(ir + k, js + l, mr)] = y * B[S(k, l, r)]; 
                }
            }
        }
    }
}  

void C_Permute(int *x, int n, int *ans) 
{
    int k = n, j;
    
    for (int i = 0; i < k; i++) {
        j = n * unif_rand();
        ans[i] = x[j];
        x[j] = x[--n];
    }
}

void C_PermuteBlock(int *x, int *table, int Ntable, int *ans)
{
    int *px, *pans;
    
    px = x;
    pans = ans;
    
    for (int j = 0; j < Ntable; j++) {
        if (table[j] > 0) {
            C_Permute(px, table[j], pans);
            px += table[j];
            pans += table[j];
        }
    }
}

void C_doPermuteBlock(int *subset, int Nsubset, int *table, int Nlevels, 
                      int *Nsubset_tmp, int *perm) 
{
    for (int i = 0; i < Nsubset; i++) Nsubset_tmp[i] = subset[i];
    C_PermuteBlock(Nsubset_tmp, table, Nlevels, perm);
}

void CR_PermuteBlockSetup(SEXP subset, SEXP block, SEXP table, SEXP orig)
{

    int Nsubset, Nlevels, *iblock, *isubset, *itable, *iorig, *itmp;
    double *subblock;

    Nsubset = LENGTH(subset);
    Nlevels = C_nlevels(block) + 1;
    C_1dtable_subset(INTEGER(block), Nlevels, INTEGER(subset), Nsubset, INTEGER(table));

    iblock = INTEGER(block);
    isubset = INTEGER(subset);
    iorig = INTEGER(orig);

    subblock = Calloc(Nsubset, double);
    itmp = Calloc(Nsubset, int);
    
    for (int i = 0; i < Nsubset; i++) {
        subblock[i] = (double) iblock[isubset[i]];
        iorig[i] = isubset[i];
    }

    rsort_with_index(subblock, iorig, Nsubset); /* first element is double */

    Free(subblock);
}

SEXP R_PermuteBlock_subset(SEXP subset, SEXP block)
{
    SEXP ans, orig, perm, table;
    int *tmp;
    
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, LENGTH(subset)));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, LENGTH(subset)));
    SET_VECTOR_ELT(ans, 2, table = allocVector(INTSXP, C_nlevels(block) + 1));
    GetRNGstate();
    
    CR_PermuteBlockSetup(subset, block, table, orig);
    tmp = Calloc(LENGTH(subset), int);

    C_doPermuteBlock(INTEGER(orig), LENGTH(subset), INTEGER(table), C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));

    Free(tmp);
                       
    PutRNGstate();
    UNPROTECT(1);
    return(ans);
}

SEXP R_PermuteBlock_weights(SEXP weights, SEXP block)
{
    SEXP ans, subset, orig, perm, table;
    int *tmp, sw, itmp;
    
    sw <- C_sum(INTEGER(weights), LENGTH(weights));
    
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 2, table = allocVector(INTSXP, C_nlevels(block) + 1));
    SET_VECTOR_ELT(ans, 3, subset = allocVector(INTSXP, sw));

    /* subset = rep(1:length(weights), weights) */
    itmp = 0;
    for (int i = 0; i < LENGTH(weights); i++) {
        for (int j = 0; j < INTEGER(weights)[i]; j++)
            INTEGER(subset)[++itmp] = i;
    }

    GetRNGstate();
    
    CR_PermuteBlockSetup(subset, block, table, orig);
    tmp = Calloc(LENGTH(subset), int);

    C_doPermuteBlock(INTEGER(orig), LENGTH(subset), INTEGER(table), C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));

    Free(tmp);
                       
    PutRNGstate();
    UNPROTECT(1);
    return(ans);
}

SEXP R_PermuteBlock_weights_subset(SEXP weights, SEXP subset, SEXP block)
{
    SEXP ans, subset2, orig, perm, table;
    int *tmp, sw, itmp;
    
    sw <- C_sum_subset(INTEGER(weights), LENGTH(weights), INTEGER(subset), LENGTH(subset));
    
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 2, table = allocVector(INTSXP, C_nlevels(block) + 1));
    SET_VECTOR_ELT(ans, 3, subset2 = allocVector(INTSXP, sw));

    /* subset = rep(subset, weights[subset]) */
    itmp = 0;
    for (int i = 0; i < LENGTH(subset); i++) {
        for (int j = 0; j < INTEGER(weights)[INTEGER(subset)[i]]; j++)
            INTEGER(subset2)[++itmp] = INTEGER(subset)[i];
    }

    GetRNGstate();
    
    CR_PermuteBlockSetup(subset2, block, table, orig);
    tmp = Calloc(LENGTH(subset2), int);

    C_doPermuteBlock(INTEGER(orig), LENGTH(subset2), INTEGER(table), C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));

    Free(tmp);
                       
    PutRNGstate();
    UNPROTECT(1);
    return(ans);
}



SEXP R_PermuteBlock(SEXP block)
{
    SEXP ans, subset, orig, perm, table;
    int *tmp, sw, itmp;

    sw = LENGTH(block);    
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, sw));
    SET_VECTOR_ELT(ans, 2, table = allocVector(INTSXP, C_nlevels(block) + 1));
    SET_VECTOR_ELT(ans, 3, subset = allocVector(INTSXP, sw));

    /* subset = 1:length(block) */
    for (int i = 0; i < sw; i++)
        INTEGER(subset)[i] = i;

    GetRNGstate();
    
    CR_PermuteBlockSetup(subset, block, table, orig);
    tmp = Calloc(LENGTH(subset), int);

    C_doPermuteBlock(INTEGER(orig), LENGTH(subset), INTEGER(table), C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));

    Free(tmp);
                       
    PutRNGstate();
    UNPROTECT(1);
    return(ans);
}

