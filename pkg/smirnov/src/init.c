
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP Smirnov_sim(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Smirnov_sim", (DL_FUNC) &Smirnov_sim, 3},
    {NULL, NULL, 0}
};

void R_init_basefun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
