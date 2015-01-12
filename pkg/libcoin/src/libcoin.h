
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "linstat.h"
#include "helpers.h"
#include "teststat.h"

SEXP R_LinstatExpCov (const SEXP data, const SEXP inputs,
                      const SEXP y, const SEXP weights);
                      