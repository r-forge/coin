
#include "libcoin.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"R_LinstatExpCov", (DL_FUNC) &R_LinstatExpCov, 4},
    {NULL, NULL, 0}
};
        
        
void R_init_libcoin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_RegisterCCallable("libcoin", "R_LinstatExpCov", (DL_FUNC) &R_LinstatExpCov);
}
