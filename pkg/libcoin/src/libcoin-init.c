
#include "libcoin.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef callMethods[] = {
    CALLDEF(R_ExpectationCovarianceStatistic, 7),
    CALLDEF(R_PermutedLinearStatistic, 8),
    CALLDEF(R_ExpectationCovarianceStatistic_2d, 9),
    CALLDEF(R_PermutedLinearStatistic_2d, 8),
    CALLDEF(R_tables, 5),
    CALLDEF(R_ChisqTest, 4),
    CALLDEF(R_MaxtypeTest, 8),
    CALLDEF(R_MaxSelectTest, 6),
    CALLDEF(R_kronecker, 2),
    {NULL, NULL, 0}
};

void R_init_libcoin
(
    DllInfo *dll
) {

    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
