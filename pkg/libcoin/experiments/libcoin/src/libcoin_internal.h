
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

#define LinearStatistic_SLOT	0
#define Expectation_SLOT	1
#define Covariance_SLOT		2
#define Variance_SLOT		2
#define ExpectationX_SLOT  	3
#define varonly_SLOT		4
#define dim_SLOT		5
#define Table_SLOT		6
