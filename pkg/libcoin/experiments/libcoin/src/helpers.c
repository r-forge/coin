
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
                  const double *B, const int r, const int s, int overwrite,
                  double *ans)
{
    int i, j, k, l, mr, js, ir;
    double y;

    if (overwrite) {
        for (i = 0; i < m * r * n * s; i++) ans[i] = 0.0;
    }

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++)
                    ans[(js + l) * mr + ir + k] += y * B[l * r + k];
            }
        }
    }
}  

void C_kronecker_sym (const double *A, const int m, 
                      const double *B, const int r, int overwrite,
                      double *ans)
{
    int i, j, k, l, mr, js, ir, s, n, tmp, mrns;
    double y;

    mr = m * r;
    n = m;
    s = r;

    if (overwrite) {
        for (i = 0; i < mr * (mr + 1) / 2; i++) ans[i] = 0.0;
    }

    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j <= i; j++) {
            js = j * s;
            y = A[S(i, j, m)]; 
            for (k = 0; k < r; k++) {
                for (l = 0; l < (j < i ? s : k + 1); l++) {
                    ans[S(ir + k, js + l, mr)] += y * B[S(k, l, r)]; 
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
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_PermuteBlock(Nsubset_tmp, table, Nlevels, perm);
}

void C_doPermute(int *subset, int Nsubset, int *Nsubset_tmp, int *perm) 
{
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_Permute(Nsubset_tmp, Nsubset, perm);
}

void C_setup_subset(int N, int *N_ans)
{
    for (int i = 0; i < N; i++) N_ans[i] = i;
}

void C_setup_subset_weights(int N, int *weights, int *sw_ans)
{
    int itmp = 0;
    /* subset = rep(1:length(weights), weights) */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < weights[i]; j++)
            sw_ans[itmp++] = i;
    }
}

void C_setup_subset_weights_subset(int Nsubset, int *weights, int *subset, int *sw_ans)
{
    /* subset = rep(subset, weights[subset]) */
    int itmp = 0;
    for (int i = 0; i < Nsubset; i++) {
        for (int j = 0; j < weights[subset[i]]; j++)
            sw_ans[++itmp] = subset[i];
    }
}

void C_order_wrt_block(int *subset, int Nsubset, int *block, int *table, int Nlevels)
{
    int *cumtable, *subset_tmp;
    
    cumtable = Calloc(Nlevels, int);
    subset_tmp = Calloc(Nsubset, int);
    Memcpy(subset_tmp, subset, Nsubset);
    cumtable[0] = 0;
    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++) cumtable[k] = cumtable[k - 1] + table[k - 1];
    
    for (int i = 0; i < Nsubset; i++)
        subset[cumtable[block[subset[i]]]++] = subset_tmp[i];
    Free(cumtable); Free(subset_tmp);
} 

SEXP R_PermuteBlock(SEXP block)
{
    SEXP ans, orig, perm;
    int N, *table, *tmp;
    
    N = LENGTH(block);    

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    table = Calloc(C_nlevels(block) + 1, int);
    C_1dtable(INTEGER(block), C_nlevels(block) + 1, LENGTH(block), table);

    C_setup_subset(N, INTEGER(orig));
    C_order_wrt_block(INTEGER(orig), N, INTEGER(block), table, C_nlevels(block) + 1);
    
    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermuteBlock(INTEGER(orig), N, table, C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));
    PutRNGstate();

    Free(tmp); Free(table);
    UNPROTECT(1);
    return(ans);
}

SEXP R_PermuteBlock_subset(SEXP subset, SEXP block)
{
    SEXP ans, orig, perm;
    int N, *table, *tmp;
    
    N = LENGTH(subset);    

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    table = Calloc(C_nlevels(block) + 1, int);
    C_1dtable_subset(INTEGER(block), C_nlevels(block) + 1, INTEGER(subset), N, table);

    Memcpy(INTEGER(orig), INTEGER(subset), N);
    C_order_wrt_block(INTEGER(orig), N, INTEGER(block), table, C_nlevels(block) + 1);
    
    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermuteBlock(INTEGER(orig), N, table, C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));
    PutRNGstate();
    Free(tmp); Free(table);
    UNPROTECT(1);
    return(ans);
}

SEXP R_PermuteBlock_weights(SEXP weights, SEXP block)
{
    SEXP ans, orig, perm;
    int N, *table, *tmp;
    
    N = C_sum(INTEGER(weights), LENGTH(weights));

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    table = Calloc(C_nlevels(block) + 1, int);
    C_1dtable_weights(INTEGER(block), C_nlevels(block) + 1, INTEGER(weights), N, table);

    C_setup_subset_weights(LENGTH(weights), INTEGER(weights), INTEGER(orig));
    C_order_wrt_block(INTEGER(orig), N, INTEGER(block), table, C_nlevels(block) + 1);
    
    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermuteBlock(INTEGER(orig), N, table, C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));
                       
    PutRNGstate();
    Free(tmp); Free(table);
    UNPROTECT(1);
    return(ans);
}

SEXP R_PermuteBlock_weights_subset(SEXP weights, SEXP subset, SEXP block)
{
    SEXP ans, orig, perm;
    int N, *table, *tmp;
    
    N = C_sum_subset(INTEGER(weights), LENGTH(weights), INTEGER(subset), LENGTH(subset));

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    table = Calloc(C_nlevels(block) + 1, int);
    C_1dtable_weights_subset(INTEGER(block), C_nlevels(block) + 1, 
                             INTEGER(weights), INTEGER(subset), 
                             LENGTH(subset), table);

    C_setup_subset_weights_subset(LENGTH(subset), INTEGER(weights), INTEGER(subset), INTEGER(orig));
    C_order_wrt_block(INTEGER(orig), N, INTEGER(block), table, C_nlevels(block) + 1);
    
    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermuteBlock(INTEGER(orig), N, table, C_nlevels(block) + 1, 
                     tmp, INTEGER(perm));
    PutRNGstate();

    Free(tmp); Free(table);
    UNPROTECT(1);
    return(ans);
}

SEXP R_Permute(SEXP n)
{
    SEXP ans, orig, perm;
    int N, *tmp;
    
    N = INTEGER(n)[0];    

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    C_setup_subset(N, INTEGER(orig));
    
    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermute(INTEGER(orig), N, tmp, INTEGER(perm));
    PutRNGstate();

    Free(tmp); 
    UNPROTECT(1);
    return(ans);
}

SEXP R_Permute_subset(SEXP subset)
{
    SEXP ans, orig, perm;
    int N, *tmp;
    
    N = LENGTH(subset);    

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    Memcpy(INTEGER(orig), INTEGER(subset), N);

    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermute(INTEGER(orig), N, tmp, INTEGER(perm));
    PutRNGstate();
    Free(tmp); 
    UNPROTECT(1);
    return(ans);
}

SEXP R_Permute_weights(SEXP weights)
{
    SEXP ans, orig, perm;
    int N, *tmp;
    
    N = C_sum(INTEGER(weights), LENGTH(weights));

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    C_setup_subset_weights(LENGTH(weights), INTEGER(weights), INTEGER(orig));

    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermute(INTEGER(orig), N, tmp, INTEGER(perm));
                       
    PutRNGstate();
    Free(tmp); 
    UNPROTECT(1);
    return(ans);
}

SEXP R_Permute_weights_subset(SEXP weights, SEXP subset)
{
    SEXP ans, orig, perm;
    int N, *tmp;
    
    N = C_sum_subset(INTEGER(weights), LENGTH(weights), INTEGER(subset), LENGTH(subset));

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, orig = allocVector(INTSXP, N));
    SET_VECTOR_ELT(ans, 1, perm = allocVector(INTSXP, N));

    C_setup_subset_weights_subset(LENGTH(subset), INTEGER(weights), INTEGER(subset), INTEGER(orig));

    tmp = Calloc(N, int);
    GetRNGstate();
    C_doPermute(INTEGER(orig), N, tmp, INTEGER(perm));
    PutRNGstate();

    Free(tmp); 
    UNPROTECT(1);
    return(ans);
}
