
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "Distributions.h"
#include "MemoryAccess.h"

SEXP R_ChisqTest(SEXP LEV, SEXP tol, SEXP give_log) 
{
    SEXP ans, MPrank, stat, pval;
    double *MPinv;
    int rank, P, Q, PQ;
    
    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;

    if (C_get_varonly(LEV))
        error("cannot compute quadratic form based on variances only");
        
    PROTECT(MPrank = R_MPinv_sym(VECTOR_ELT(LEV, Covariance_SLOT), tol));
    MPinv = REAL(VECTOR_ELT(MPrank, 0));
    rank = INTEGER(VECTOR_ELT(MPrank, 1))[0];
        
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    
    REAL(stat)[0] = C_quadform(PQ, C_get_LinearStatistic(LEV),
                               C_get_Expectation(LEV), MPinv);
    REAL(pval)[0] = C_chisq_pvalue(REAL(stat)[0], rank,
                                   INTEGER(give_log)[0]);
    UNPROTECT(2);
    return(ans);
}

SEXP R_MaxabsstatTest(SEXP LEV, SEXP tol, SEXP give_log, SEXP maxpts, SEXP releps, SEXP abseps)
{
    SEXP ans, stat, pval;
    int P, Q, PQ;

    P = C_get_P(LEV);
    Q = C_get_Q(LEV);
    PQ = P * Q;
            
    if (C_get_varonly(LEV) && PQ > 1)
            error("cannot compute adjusted p-value based on variances only");
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));

    REAL(stat)[0] =  C_maxabsstat_Covariance(PQ, C_get_LinearStatistic(LEV), 
                                             C_get_Expectation(LEV), 
                                             C_get_Covariance(LEV), REAL(tol)[0]);
    REAL(pval)[0] = C_maxabsstat_pvalue(REAL(stat)[0], C_get_Covariance(LEV),
                                        PQ,
                                        INTEGER(maxpts)[0], REAL(releps)[0], 
                                        REAL(abseps)[0], REAL(tol)[0]);
                                        
    if (INTEGER(give_log)[0]) REAL(pval)[0] = log(REAL(pval)[0]);
    
    UNPROTECT(1);
    return(ans);
}
