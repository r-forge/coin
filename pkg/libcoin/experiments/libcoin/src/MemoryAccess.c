
#include "libcoin_internal.h"

int C_get_P(SEXP LECV) 
{
    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[0]);
}

int C_get_Q(SEXP LECV) 
{
    return(INTEGER(VECTOR_ELT(LECV, dim_SLOT))[1]);
}

int C_get_varonly(SEXP LECV)
{
    return(INTEGER(VECTOR_ELT(LECV, varonly_SLOT))[0]);
}

double* C_get_LinearStatistic(SEXP LECV)
{
    REAL(VECTOR_ELT(LECV, LinearStatistic_SLOT));
}

double* C_get_Expectation(SEXP LECV)
{
    REAL(VECTOR_ELT(LECV, Expectation_SLOT));
}

double* C_get_Covariance(SEXP LECV)
{
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract covariance from variance only object");
    REAL(VECTOR_ELT(LECV, Covariance_SLOT));
}

double* C_get_Variance(SEXP LECV)
{
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (!C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract variance from covariance object");
    REAL(VECTOR_ELT(LECV, Variance_SLOT));
}

double* C_get_ExpectationX(SEXP LECV)
{
    REAL(VECTOR_ELT(LECV, ExpectationX_SLOT));
}

int* C_get_Table(SEXP LECV)
{
    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    INTEGER(VECTOR_ELT(LECV, Table_SLOT));
}

int* C_get_dimTable (SEXP LECV) {
    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");              
    INTEGER(getAttrib(VECTOR_ELT(LECV, Table_SLOT), R_DimSymbol));
}

SEXP R_init_LECV(SEXP P, SEXP Q, SEXP varonly)
{
    SEXP ans, vo, d, names;
    int p, q, pq;
    
    if (!isInteger(P) || LENGTH(P) != 1)
        error("P is not a scalar integer");
    if (INTEGER(P)[0] <= 0)
        error("P is not positive");

    if (!isInteger(Q) || LENGTH(Q) != 1)
        error("Q is not a scalar integer");
    if (INTEGER(Q)[0] <= 0)
        error("Q is not positive");

    if (!isInteger(varonly) || LENGTH(varonly) != 1)
        error("varonly is not a scalar integer");
    if (INTEGER(varonly)[0] < 0 || INTEGER(varonly)[0] > 1)
            error("varonly is not 0 or 1");

    p = INTEGER(P)[0];
    q = INTEGER(Q)[0];
    pq = p * q;

    PROTECT(ans = allocVector(VECSXP, dim_SLOT + 1));
    PROTECT(names = allocVector(STRSXP, dim_SLOT + 1)); 

    SET_VECTOR_ELT(ans, LinearStatistic_SLOT, 
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, LinearStatistic_SLOT, 
                   mkChar("LinearStatistic"));

    SET_VECTOR_ELT(ans, Expectation_SLOT, 
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, Expectation_SLOT, 
                   mkChar("Expectation"));
                   
    SET_VECTOR_ELT(ans, varonly_SLOT, 
                   vo = allocVector(INTSXP, 1));
    SET_STRING_ELT(names, varonly_SLOT, 
                   mkChar("varonly"));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, Variance_SLOT, 
                       allocVector(REALSXP, pq));
        SET_STRING_ELT(names, Variance_SLOT, 
                       mkChar("Variance"));
        INTEGER(vo)[0] = 1;
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, Covariance_SLOT, 
                       mkChar("Covariance"));
        INTEGER(vo)[0] = 0;
    }

    SET_VECTOR_ELT(ans, ExpectationX_SLOT, 
                   allocVector(REALSXP, p));
    SET_STRING_ELT(names, ExpectationX_SLOT,
                   mkChar("ExpectationX"));

    SET_VECTOR_ELT(ans, dim_SLOT, 
                   d = allocVector(INTSXP, 2));
    SET_STRING_ELT(names, dim_SLOT, 
                   mkChar("dimension"));
    INTEGER(d)[0] = p;
    INTEGER(d)[1] = q;
    
    namesgets(ans, names);
    
    UNPROTECT(2);
    return(ans);
}

SEXP R_init_LECV_2d(SEXP P, SEXP Q, SEXP varonly, SEXP Lx, SEXP Ly, SEXP Lb) 
{
    SEXP ans, vo, d, names, tab, tabdim;
    int p, q, pq;
    
    if (!isInteger(P) || LENGTH(P) != 1)
        error("P is not a scalar integer");
    if (INTEGER(P)[0] <= 0)
        error("P is not positive");

    if (!isInteger(Q) || LENGTH(Q) != 1)
        error("Q is not a scalar integer");
    if (INTEGER(Q)[0] <= 0)
        error("Q is not positive");

    if (!isInteger(varonly) || LENGTH(varonly) != 1)
        error("varonly is not a scalar integer");
    if (INTEGER(varonly)[0] < 0 || INTEGER(varonly)[0] > 1)
            error("varonly is not 0 or 1");

    if (!isInteger(Lx) || LENGTH(Lx) != 1)
        error("Lx is not a scalar integer");
    if (INTEGER(Lx)[0] <= 0)
        error("Lx is not positive");

    if (!isInteger(Ly) || LENGTH(Ly) != 1)
        error("Ly is not a scalar integer");
    if (INTEGER(Ly)[0] <= 0)
        error("Ly is not positive");

    if (!isInteger(Lb) || LENGTH(Lb) != 1)
        error("Lb is not a scalar integer");
    if (INTEGER(Lb)[0] <= 0)
        error("Lb is not positive");

    p = INTEGER(P)[0];
    q = INTEGER(Q)[0];
    pq = p * q;

    PROTECT(ans = allocVector(VECSXP, Table_SLOT + 1));
    PROTECT(names = allocVector(STRSXP, Table_SLOT + 1)); 

    SET_VECTOR_ELT(ans, LinearStatistic_SLOT, 
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, LinearStatistic_SLOT, 
                   mkChar("LinearStatistic"));

    SET_VECTOR_ELT(ans, Expectation_SLOT, 
                   allocVector(REALSXP, pq));
    SET_STRING_ELT(names, Expectation_SLOT, 
                   mkChar("Expectation"));
                   
    SET_VECTOR_ELT(ans, varonly_SLOT, 
                   vo = allocVector(INTSXP, 1));
    SET_STRING_ELT(names, varonly_SLOT, 
                   mkChar("varonly"));
    if (INTEGER(varonly)[0]) {
        SET_VECTOR_ELT(ans, Variance_SLOT, 
                       allocVector(REALSXP, pq));
        SET_STRING_ELT(names, Variance_SLOT, 
                       mkChar("Variance"));
        INTEGER(vo)[0] = 1;
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, Covariance_SLOT, 
                       mkChar("Covariance"));
        INTEGER(vo)[0] = 0;
    }

    SET_VECTOR_ELT(ans, ExpectationX_SLOT, 
                   allocVector(REALSXP, p));
    SET_STRING_ELT(names, ExpectationX_SLOT,
                   mkChar("ExpectationX"));
                                          
    SET_VECTOR_ELT(ans, dim_SLOT, 
                   d = allocVector(INTSXP, 2));
    SET_STRING_ELT(names, dim_SLOT, 
                   mkChar("dimension"));
    INTEGER(d)[0] = p;
    INTEGER(d)[1] = q;

    PROTECT(tabdim = allocVector(INTSXP, 3));                   
    INTEGER(tabdim)[0] = INTEGER(Lx)[0] + 1;
    INTEGER(tabdim)[1] = INTEGER(Ly)[0] + 1;
    INTEGER(tabdim)[2] = INTEGER(Lb)[0];
    SET_VECTOR_ELT(ans, Table_SLOT, 
                   tab = allocVector(INTSXP, 
                       INTEGER(tabdim)[0] * 
                       INTEGER(tabdim)[1] *
                       INTEGER(tabdim)[2]));
    dimgets(tab, tabdim);
    SET_STRING_ELT(names, Table_SLOT, 
                   mkChar("Table"));

    namesgets(ans, names);
    
    UNPROTECT(3);
    return(ans);
}

