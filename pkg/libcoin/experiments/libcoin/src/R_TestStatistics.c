
#include "libcoin.h"
#include "TestStatistics.h"

SEXP R_quadform(SEXP LinearStatistic, SEXP Expectation, SEXP MPinv)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = CR_quadform(LinearStatistic, Expectation, MPinv);
    UNPROTECT(1);
    return(ans);
}
