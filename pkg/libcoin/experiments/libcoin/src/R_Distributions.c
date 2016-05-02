
#include "libcoin.h"
#include "Tables.h"
#include "Sums.h"
#include "Distributions.h"
#include "helpers.h"

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
    C_1dtable_weights(INTEGER(block), C_nlevels(block) + 1, INTEGER(weights), 
                      LENGTH(weights), table);

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
