
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

/* StreitbergRoehmel.c */
extern SEXP C_StreitbergRoehmel
(
    SEXP score,
    SEXP m
);

static const R_CallMethodDef callMethods[] = {
    CALLDEF(C_StreitbergRoehmel, 2),
    {NULL, NULL, 0}
};

void attribute_visible R_init_coin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
