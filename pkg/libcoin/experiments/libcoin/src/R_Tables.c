
#include "libcoin.h"
#include "Tables.h"
#include "helpers.h"

void CR_2dtable(SEXP x, SEXP y, int *ans)
{
    C_2dtable(INTEGER(x), C_nlevels(x) + 1,
              INTEGER(y), C_nlevels(y) + 1,
              LENGTH(x), ans);
}

SEXP R_2dtable(SEXP x, SEXP y)
{
    SEXP ans;

    PROTECT(ans = allocMatrix(INTSXP, C_nlevels(x) + 1, C_nlevels(y) + 1)); 
    CR_2dtable(x, y, INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_2dtable_subset(SEXP x, SEXP y, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(INTSXP, C_nlevels(x) + 1, C_nlevels(y) + 1)); 
    C_2dtable_subset(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                     INTEGER(subset), LENGTH(subset), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_2dtable_weights(SEXP x, SEXP y, SEXP weights) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(INTSXP, C_nlevels(x) + 1, C_nlevels(y) + 1)); 
    C_2dtable_weights(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                      INTEGER(weights), LENGTH(x), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_2dtable_weights_subset(SEXP x, SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans;
    
    PROTECT(ans = allocMatrix(INTSXP, C_nlevels(x) + 1, C_nlevels(y) + 1)); 
    C_2dtable_weights_subset(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                             INTEGER(weights), INTEGER(subset), LENGTH(subset), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}


SEXP R_2dtable_block(SEXP x, SEXP y, SEXP block)
{
    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, (C_nlevels(x) + 1) * (C_nlevels(y) + 1) * C_nlevels(block))); 
    C_2dtable_block(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                    INTEGER(block), C_nlevels(block), LENGTH(x), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}


SEXP R_2dtable_subset_block(SEXP x, SEXP y, SEXP subset, SEXP block) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(INTSXP, (C_nlevels(x) + 1) * (C_nlevels(y) + 1) * C_nlevels(block))); 
    C_2dtable_subset_block(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                     INTEGER(subset), LENGTH(subset), INTEGER(block), C_nlevels(block), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_2dtable_weights_block(SEXP x, SEXP y, SEXP weights, SEXP block) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(INTSXP, (C_nlevels(x) + 1) * (C_nlevels(y) + 1) * C_nlevels(block))); 
    C_2dtable_weights_block(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                      INTEGER(weights), INTEGER(block), C_nlevels(block), LENGTH(x), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_2dtable_weights_subset_block(SEXP x, SEXP y, SEXP weights, SEXP subset, SEXP block) 
{
    SEXP ans;
    
    PROTECT(ans = allocVector(INTSXP, (C_nlevels(x) + 1) * (C_nlevels(y) + 1) * C_nlevels(block))); 
    C_2dtable_weights_subset_block(INTEGER(x), C_nlevels(x) + 1, INTEGER(y), C_nlevels(y) + 1, 
                             INTEGER(weights), INTEGER(subset), LENGTH(subset), INTEGER(block), C_nlevels(block), 
                             INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}


SEXP R_1dtable(SEXP y) 
{
    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, C_nlevels(y) + 1)); 
    C_1dtable(INTEGER(y), C_nlevels(y) + 1, LENGTH(y), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_1dtable_subset(SEXP y, SEXP subset) 
{
    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, C_nlevels(y) + 1)); 
    C_1dtable_subset(INTEGER(y), C_nlevels(y) + 1, INTEGER(subset), 
                     LENGTH(subset), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_1dtable_weights(SEXP y, SEXP weights) 
{
    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, C_nlevels(y) + 1)); 
    C_1dtable_subset(INTEGER(y), C_nlevels(y) + 1, INTEGER(weights), 
                     LENGTH(weights), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

SEXP R_1dtable_weights_subset(SEXP y, SEXP weights, SEXP subset) 
{
    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, C_nlevels(y) + 1)); 
    C_1dtable_weights_subset(INTEGER(y), C_nlevels(y) + 1, 
                             INTEGER(weights), INTEGER(subset), 
                             LENGTH(subset), INTEGER(ans));
    UNPROTECT(1);
    return(ans);
}

/*

SEXP R_d2s(SEXP x) {

    SEXP ans;
    int nr, nc, i, j, s, *ix, cnt, *ians;
    
    nr = NROW(x);
    nc = NCOL(x);
    ix = INTEGER(x),
    
    s = 0;
    for (i = 0; i < nr * nc; i++)
        if (ix[i] != 0) s++;
        
//    if (s / (nr + nc) > .5)
    PROTECT(ans = allocMatrix(INTSXP, s, 3));
    ians = INTEGER(ans);
    
    cnt = 0;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            if (ix[i + j * nr] != 0) {
                ians[cnt] = i;
                ians[NROW(ans) + cnt] = j;
                ians[2 * NROW(ans) + cnt] = ix[i + j * nr];
                cnt++;
            }
       }
     }
     UNPROTECT(1);
     return(ans);
}

*/