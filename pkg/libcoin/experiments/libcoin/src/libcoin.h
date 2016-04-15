
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */

#define ALTERNATIVE_twosided               1    
#define ALTERNATIVE_less                   2    
#define ALTERNATIVE_greater                3    
