
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

SEXP R_maxtype(SEXP LinearStatistic, SEXP Expectation, SEXP CoVariance, SEXP type, SEXP tol)
{
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    if (INTEGER(type)[0] == ALTERNATIVE_twosided)
        REAL(ans)[0] = CR_maxabsstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_less)
        REAL(ans)[0] = CR_maxstat(LinearStatistic, Expectation, CoVariance, tol);
    if (INTEGER(type)[0] == ALTERNATIVE_greater)
        REAL(ans)[0] = CR_minstat(LinearStatistic, Expectation, CoVariance, tol);
    UNPROTECT(1);
    return(ans);
}
