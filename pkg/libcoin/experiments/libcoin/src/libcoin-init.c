
#include "libcoin.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"R_ExpectationCovarianceStatistic", (DL_FUNC) &R_ExpectationCovarianceStatistic, 6},
    {"R_PermutedLinearStatistic", (DL_FUNC) &R_PermutedLinearStatistic, 7},
    {"R_ExpectationCovarianceStatistic_2d", (DL_FUNC) &R_ExpectationCovarianceStatistic_2d, 8},
    {"R_PermutedLinearStatistic_2d", (DL_FUNC) &R_PermutedLinearStatistic_2d, 7},
    {NULL, NULL, 0}
};
        
void R_init_libcoin(DllInfo *dll)
{ 
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
    R_RegisterCCallable("mvtnorm", "R_ExpectationCovarianceStatistic", (DL_FUNC) &R_ExpectationCovarianceStatistic);
    R_RegisterCCallable("mvtnorm", "R_PermutedLinearStatistic", (DL_FUNC) &R_PermutedLinearStatistic);
    R_RegisterCCallable("mvtnorm", "R_ExpectationCovarianceStatistic_2d", (DL_FUNC) &R_ExpectationCovarianceStatistic_2d);
    R_RegisterCCallable("mvtnorm", "R_PermutedLinearStatistic_2d", (DL_FUNC) &R_PermutedLinearStatistic_2d);
}
