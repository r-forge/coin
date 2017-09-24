
#include <Rinternals.h>

extern SEXP R_ExpectationCovarianceStatistic
(
    const SEXP x,
    const SEXP y,
    const SEXP weights,
    const SEXP subset,
    const SEXP block,
    const SEXP varonly,
    const SEXP tol
);
