
#include "libcoin_internal.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "Distributions.h"

SEXP R_ChisqTest(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv, 
                 SEXP rank, SEXP give_log)
{
    SEXP ans, stat, pval;
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    
    REAL(stat)[0] = C_quadform(LENGTH(LinearStatistic), REAL(LinearStatistic), 
                               REAL(Expectation), REAL(MPinv));
    REAL(pval)[0] = C_chisq_pvalue(REAL(stat)[0], INTEGER(rank)[0],
                                   INTEGER(give_log)[0]);
    UNPROTECT(1);
    return(ans);
}

SEXP R_MaxabsstatTest(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, 
                      SEXP give_log, SEXP tol, SEXP maxpts, SEXP releps, SEXP abseps)
{
    SEXP ans, stat, pval;
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, stat = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, pval = allocVector(REALSXP, 1));
    
    REAL(stat)[0] =  CR_maxabsstat(LinearStatistic, Expectation, CoVariance, tol);
    REAL(pval)[0] = C_maxabsstat_pvalue(REAL(stat)[0], REAL(CoVariance), 
                                        LENGTH(LinearStatistic),
                                        INTEGER(maxpts)[0], INTEGER(releps)[0], 
                                        INTEGER(abseps)[0], REAL(tol)[0]);
    if (INTEGER(give_log)[0]) REAL(pval)[0] = log(REAL(pval)[0]);
    UNPROTECT(1);
    return(ans);
}

