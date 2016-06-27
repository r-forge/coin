
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

#define LinearStatistic_SLOT		0	
#define Expectation_SLOT		1
#define Covariance_SLOT			2
#define Variance_SLOT			2
#define MPinv_SLOT			3
#define ExpectationX_SLOT  		4
#define varonly_SLOT			5
#define dim_SLOT			6
#define ExpectationInfluence_SLOT	7
#define CovarianceInfluence_SLOT	8
#define Work_SLOT			9
#define TableBlock_SLOT			10
#define Sumweights_SLOT			11
#define Table_SLOT			12

SEXP R_MPinv_sym (SEXP x, SEXP tol);