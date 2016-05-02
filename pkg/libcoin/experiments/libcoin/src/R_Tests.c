
#include "libcoin.h"
#include "Utils.h"
#include "TestStatistics.h"
#include "Distributions.h"
#include "R_Utils.h"

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
