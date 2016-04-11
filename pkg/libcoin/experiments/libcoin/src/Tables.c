
#include <R.h>
#include <Rinternals.h>

void C_2int_set_zero(int *ians, int inx, int iny) {

    int j, k;
    
    for (j = 0; j < inx; j++) {
         for (k = 0; k < iny; k++) {
              ians[k * inx + j] = 0;
         }
    }
}

void C_2int_table(int *ix, int inx, int *iy, int iny, int n, int *ians) {

    C_2int_set_zero(ians, inx, iny);
    
    for (int i = 0; i < n; i++)  
        ians[ix[i] + iy[i] * inx]++;

}

SEXP R_2int_table(SEXP x, SEXP nx, SEXP y, SEXP ny) {

    SEXP ans;

    PROTECT(ans = allocMatrix(INTSXP, INTEGER(nx)[0], INTEGER(ny)[0])); 
    
    C_2int_table(INTEGER(x), INTEGER(nx)[0], 
                 INTEGER(y), INTEGER(ny)[0], 
                 LENGTH(x), INTEGER(ans));
    
    UNPROTECT(1);
    return(ans);
}


SEXP R_2int_s_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP subset) {

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
    
    C_2int_set_zero(ians, inx, iny);
    
    for (i = 0; i < n; i++)  
        ians[ix[is[i]] + iy[is[i]] * inx]++;

    UNPROTECT(1);
    return(ans);
}


SEXP R_2int_w_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP weights) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy, *iw, idx;
    
    n = LENGTH(x); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    iw = INTEGER(weights);
    
    C_2int_set_zero(ians, inx, iny);
    
    for (i = 0; i < n; i++) {
        idx = ix[i] + iy[i] * inx;
        ians[idx] = ians[idx] + iw[i];
    }

    UNPROTECT(1);
    return(ans);
}

SEXP R_2int_s_w_table(SEXP x, SEXP nx, SEXP y, SEXP ny, SEXP subset, SEXP weights) {

    SEXP ans;
    int n, inx, iny, i, *ians, *ix, *iy, *is, *iw, idx;
    
    n = LENGTH(subset); 
    inx = INTEGER(nx)[0];
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocMatrix(INTSXP, inx, iny)); 
    ians = INTEGER(ans);
    ix = INTEGER(x);
    iy = INTEGER(y);
    is = INTEGER(subset);
    iw = INTEGER(weights);
    
    C_2int_set_zero(ians, inx, iny);
   
    for (i = 0; i < n; i++) {
        idx = ix[is[i]] + iy[is[i]] * inx;
        ians[idx] = ians[idx] + iw[is[i]];
    }
        
    UNPROTECT(1);
    return(ans);
}


void C_1int_set_zero(int *ians, int iny) {

    int k;
    
     for (k = 0; k < iny; k++)
          ians[k] = 0;
}

void C_1int_table(int *iy, int iny, int n, int *ians) {

    C_1int_set_zero(ians, iny);
    
    for (int i = 0; i < n; i++)  
        ians[iy[i]]++;

}

SEXP R_1int_table(SEXP y, SEXP ny) {

    SEXP ans;

    PROTECT(ans = allocVector(INTSXP, INTEGER(ny)[0])); 
    
    C_1int_table(INTEGER(y), INTEGER(ny)[0], 
                 LENGTH(y), INTEGER(ans));
    
    UNPROTECT(1);
    return(ans);
}


SEXP R_1int_s_table(SEXP y, SEXP ny, SEXP subset) {

    SEXP ans;
    int n, iny, i, *ians, *iy, *is;
    
    n = LENGTH(subset); 
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocVector(INTSXP, iny)); 
    ians = INTEGER(ans);
    iy = INTEGER(y);
    is = INTEGER(subset);
    
    C_1int_set_zero(ians, iny);
    
    for (i = 0; i < n; i++)  
        ians[iy[is[i]]]++;

    UNPROTECT(1);
    return(ans);
}


SEXP R_1int_w_table(SEXP y, SEXP ny, SEXP weights) {

    SEXP ans;
    int n, iny, i, *ians, *iy, *iw, idx;
    
    n = LENGTH(y); 
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocVector(INTSXP, iny)); 
    ians = INTEGER(ans);
    iy = INTEGER(y);
    iw = INTEGER(weights);
    
    C_1int_set_zero(ians, iny);
    
    for (i = 0; i < n; i++) {
        idx = iy[i];
        ians[idx] = ians[idx] + iw[i];
    }

    UNPROTECT(1);
    return(ans);
}

SEXP R_1int_s_w_table(SEXP y, SEXP ny, SEXP subset, SEXP weights) {

    SEXP ans;
    int n, iny, i, *ians, *iy, *is, *iw, idx;
    
    n = LENGTH(subset); 
    iny = INTEGER(ny)[0];
    
    PROTECT(ans = allocVector(INTSXP, iny)); 
    ians = INTEGER(ans);
    iy = INTEGER(y);
    is = INTEGER(subset);
    iw = INTEGER(weights);
    
    C_1int_set_zero(ians, iny);
   
    for (i = 0; i < n; i++) {
        idx = iy[is[i]];
        ians[idx] = ians[idx] + iw[is[i]];
    }
        
    UNPROTECT(1);
    return(ans);
}

int nrow(SEXP x) {
    SEXP a;
    
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(LENGTH(x));
    } else {
        return(INTEGER(getAttrib(x, R_DimSymbol))[0]);
    }
}

int ncol(SEXP x) {
    SEXP a;
    
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) {
        return(LENGTH(x));
    } else {
        return(INTEGER(getAttrib(x, R_DimSymbol))[1]);
    }
}


SEXP R_d2s(SEXP x) {

    SEXP ans;
    int nr, nc, i, j, s, *ix, cnt, *ians;
    
    nr = nrow(x);
    nc = ncol(x);
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
                ians[nrow(ans) + cnt] = j;
                ians[2 * nrow(ans) + cnt] = ix[i + j * nr];
                cnt++;
            }
       }
     }
     UNPROTECT(1);
     return(ans);
}
