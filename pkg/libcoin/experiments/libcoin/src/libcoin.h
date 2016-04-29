
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h> /* for dgesdd */

#define ALTERNATIVE_twosided               1    
#define ALTERNATIVE_less                   2    
#define ALTERNATIVE_greater                3    

/* S[i, j] for n x n symmetric matrix in lower packed storage allowing for i < j */
#define S(i, j, n) ((i) >= (j) ? (n) * (j) + (i) - (j) * ((j) + 1) / 2 : (n) * (i) + (j) - (i) * ((i) + 1) / 2)
