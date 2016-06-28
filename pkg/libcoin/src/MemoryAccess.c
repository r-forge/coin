
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
    return(REAL(VECTOR_ELT(LECV, LinearStatistic_SLOT)));
}

double* C_get_Expectation(SEXP LECV)
{
    return(REAL(VECTOR_ELT(LECV, Expectation_SLOT)));
}

double* C_get_Covariance(SEXP LECV)
{
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract covariance from variance only object");
    return(REAL(VECTOR_ELT(LECV, Covariance_SLOT)));
}

double* C_get_MPinv(SEXP LECV)
{
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract MPinv from variance only object");
    return(REAL(VECTOR_ELT(LECV, MPinv_SLOT)));
}

double* C_get_Variance(SEXP LECV)
{
    int PQ = C_get_P(LECV) * C_get_Q(LECV);
    if (!C_get_varonly(LECV) && PQ > 1)
        error("Cannot extract variance from covariance object");
    return(REAL(VECTOR_ELT(LECV, Variance_SLOT)));
}

double* C_get_ExpectationX(SEXP LECV)
{
    return(REAL(VECTOR_ELT(LECV, ExpectationX_SLOT)));
}

double* C_get_ExpectationInfluence(SEXP LECV)
{
    return(REAL(VECTOR_ELT(LECV, ExpectationInfluence_SLOT)));
}

double* C_get_CovarianceInfluence(SEXP LECV)
{
    return(REAL(VECTOR_ELT(LECV, CovarianceInfluence_SLOT)));
}

double* C_get_Work(SEXP LECV)
{
    return(REAL(VECTOR_ELT(LECV, Work_SLOT)));
}

int* C_get_TableBlock(SEXP LECV)
{
    return(INTEGER(VECTOR_ELT(LECV, TableBlock_SLOT)));
}

int* C_get_Sumweights(SEXP LECV)
{
    return(INTEGER(VECTOR_ELT(LECV, Sumweights_SLOT)));
}

int* C_get_Table(SEXP LECV)
{
    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");
    return(INTEGER(VECTOR_ELT(LECV, Table_SLOT)));
}

int* C_get_dimTable (SEXP LECV) {
    if (LENGTH(LECV) <= Table_SLOT)
        error("Cannot extract table from object");              
    return(INTEGER(getAttrib(VECTOR_ELT(LECV, Table_SLOT), 
                             R_DimSymbol)));
}

SEXP R_init_LECV(SEXP P, SEXP Q, SEXP varonly, SEXP Lb)
{
    SEXP ans, vo, d, names;
    int p, q, pq, lb;
    
    if (!isInteger(P) || LENGTH(P) != 1)
        error("P is not a scalar integer");
    if (INTEGER(P)[0] <= 0)
        error("P is not positive");

    if (!isInteger(Q) || LENGTH(Q) != 1)
        error("Q is not a scalar integer");
    if (INTEGER(Q)[0] <= 0)
        error("Q is not positive");

    if (!isInteger(Lb) || LENGTH(Lb) != 1)
        error("Lb is not a scalar integer");
    if (INTEGER(Lb)[0] <= 0)
        error("Lb is not positive");

    if (!isInteger(varonly) || LENGTH(varonly) != 1)
        error("varonly is not a scalar integer");
    if (INTEGER(varonly)[0] < 0 || INTEGER(varonly)[0] > 1)
            error("varonly is not 0 or 1");

    p = INTEGER(P)[0];
    q = INTEGER(Q)[0];
    lb = INTEGER(Lb)[0];
    pq = p * q;

    /* Table_SLOT is always last and only used in 2d case 
       ie omitted here */
    PROTECT(ans = allocVector(VECSXP, Table_SLOT));
    PROTECT(names = allocVector(STRSXP, Table_SLOT)); 

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
        
        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, lb * q));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT, 
                       mkChar("VarianceInfluence"));
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 3 * p + 1));
        SET_STRING_ELT(names, Work_SLOT, 
                       mkChar("Work"));
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, Covariance_SLOT, 
                       mkChar("Covariance"));
        INTEGER(vo)[0] = 0;
        SET_VECTOR_ELT(ans, MPinv_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, MPinv_SLOT, 
                       mkChar("MPinv"));

        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, lb * q * (q + 1)));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT, 
                       mkChar("CovarianceInfluence"));
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 
                           p + 2 * p * (p + 1) / 2 + 1));
        SET_STRING_ELT(names, Work_SLOT, 
                       mkChar("Work"));
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
    
    SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                   allocVector(REALSXP, lb * q));
    SET_STRING_ELT(names, ExpectationInfluence_SLOT, 
                   mkChar("ExpectationInfluence"));
                   
    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(INTSXP, lb + 1));
    SET_STRING_ELT(names, TableBlock_SLOT, 
                   mkChar("TableBlock"));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(INTSXP, lb));
    SET_STRING_ELT(names, Sumweights_SLOT, 
                   mkChar("Sumweights"));
    
    namesgets(ans, names);
    
    UNPROTECT(2);
    return(ans);
}

SEXP R_init_LECV_2d(SEXP P, SEXP Q, SEXP varonly, SEXP Lx, SEXP Ly, SEXP Lb) 
{
    SEXP ans, vo, d, names, tab, tabdim;
    int p, q, pq, lb;
    
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
    lb = INTEGER(Lb)[0];
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

        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, lb * q));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT, 
                       mkChar("VarianceInfluence"));
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 2 * p));
        SET_STRING_ELT(names, Work_SLOT, 
                       mkChar("Work"));
    } else  {
        SET_VECTOR_ELT(ans, Covariance_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, Covariance_SLOT, 
                       mkChar("Covariance"));
        INTEGER(vo)[0] = 0;
        SET_VECTOR_ELT(ans, MPinv_SLOT, 
                       allocVector(REALSXP, 
                                   pq * (pq + 1) / 2));
        SET_STRING_ELT(names, MPinv_SLOT, 
                       mkChar("MPinv"));

        SET_VECTOR_ELT(ans, CovarianceInfluence_SLOT,
                       allocVector(REALSXP, lb * q * (q + 1)));
        SET_STRING_ELT(names, CovarianceInfluence_SLOT, 
                       mkChar("CovarianceInfluence"));
        SET_VECTOR_ELT(ans, Work_SLOT,
                       allocVector(REALSXP, 
                           2 * p * (p + 1) / 2 + 1));
        SET_STRING_ELT(names, Work_SLOT, 
                       mkChar("Work"));
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
    
    SET_VECTOR_ELT(ans, ExpectationInfluence_SLOT,
                   allocVector(REALSXP, lb * q));
    SET_STRING_ELT(names, ExpectationInfluence_SLOT, 
                   mkChar("ExpectationInfluence"));
                   
    SET_VECTOR_ELT(ans, TableBlock_SLOT,
                   allocVector(INTSXP, lb + 1));
    SET_STRING_ELT(names, TableBlock_SLOT, 
                   mkChar("TableBlock"));

    SET_VECTOR_ELT(ans, Sumweights_SLOT,
                   allocVector(INTSXP, lb));
    SET_STRING_ELT(names, Sumweights_SLOT, 
                   mkChar("Sumweights"));

    PROTECT(tabdim = allocVector(INTSXP, 3));                   
    INTEGER(tabdim)[0] = INTEGER(Lx)[0] + 1;
    INTEGER(tabdim)[1] = INTEGER(Ly)[0] + 1;
    INTEGER(tabdim)[2] = lb;
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

