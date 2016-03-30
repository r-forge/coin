
#include <R.h>
#include <Rinternals.h>

void set_zero(int *ians, int inx, int iny) {

    int j, k;
    
    for (j = 0; j < inx; j++) {
         for (k = 0; k < iny; k++) {
              ians[k * inx + j] = 0;
         }
    }
}

SEXP R_int_table(SEXP x, SEXP nx, SEXP y, SEXP ny) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy;
    
    n = LENGTH(x); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    
    set_zero(ians, inx, iny);
    
    for (i = 0; i < n; i++)  
        ians[ix[i] + iy[i] * inx]++;

    UNPROTECT(1);
    return(ans);
}

SEXP R_int_s_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP subset) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy, *is;
    
    n = LENGTH(subset); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    is = INTEGER(subset);
    
    set_zero(ians, inx, iny);
    
    for (i = 0; i < n; i++)  
        ians[ix[is[i]] + iy[is[i]] * inx]++;

    UNPROTECT(1);
    return(ans);
}


SEXP R_int_w_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP weights) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy, *iw;
    
    n = LENGTH(x); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    iw = INTEGER(weights);
    
    set_zero(ians, inx, iny);
    
    for (i = 0; i < n; i++)  
        ians[ix[i] + iy[i] * inx] = ians[ix[i] + iy[i] * inx] + iw[i];

    UNPROTECT(1);
    return(ans);
}

SEXP R_int_s_w_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP subset, SEXP weights) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy, *is, *iw;
    
    n = LENGTH(subset); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    is = INTEGER(subset);
    iw = INTEGER(weights);
    
    set_zero(ians, inx, iny);
   
    for (i = 0; i < n; i++)  
        ians[ix[is[i]] + iy[is[i]] * inx] = ians[ix[is[i]] + iy[is[i]] * inx] + iw[is[i]];
        
    UNPROTECT(1);
    return(ans);
}
