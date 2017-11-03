
/* C Header */

/*
    TO NOT EDIT THIS FILE
    
    Edit `libcoin.w' and run `nuweb -r libcoin.w'
*/

#include "libcoin_internal.h"
/* Function Prototypes */

extern 
       SEXP R_ExpectationCovarianceStatistic
       (
       /* User Interface Inputs */

       /* R x Input */

           SEXP x,
       
       /* R y Input */

           SEXP y,
       
       /* R weights Input */

           SEXP weights
       ,
       /* R subset Input */

           SEXP subset
       ,
       /* R block Input */

           SEXP block
       ,
       
       SEXP varonly,
       SEXP tol
       )
       ;
extern 
       SEXP R_PermutedLinearStatistic
       (
           /* User Interface Inputs */
           
           /* R x Input */

               SEXP x,
           
           /* R y Input */

               SEXP y,
           
           /* R weights Input */

               SEXP weights
           ,
           /* R subset Input */

               SEXP subset
           ,
           /* R block Input */

               SEXP block
           ,
           
           SEXP nperm
       )
       ;
extern 
       SEXP R_StandardisePermutedLinearStatistic
       (
           SEXP LECV
       )
       ;
extern 
       SEXP R_ExpectationCovarianceStatistic_2d
       (
       /* 2d User Interface Inputs */

       /* R x Input */

           SEXP x,
       
       SEXP ix,
       /* R y Input */

           SEXP y,
       
       SEXP iy,
       /* R weights Input */

           SEXP weights
       ,
       /* R subset Input */

           SEXP subset
       ,
       /* R block Input */

           SEXP block
       ,
       
       SEXP varonly,
       SEXP tol
       )
       ;
extern 
       SEXP R_PermutedLinearStatistic_2d
       (
           /* R x Input */
           
               SEXP x,
           
           SEXP ix,
           /* R y Input */
           
               SEXP y,
           
           SEXP iy,
           /* R block Input */
           
               SEXP block
           ,
           SEXP nperm,
           SEXP itable
       )
       ;
extern 
       SEXP R_QuadraticTest
       (
           /* R LECV Input */
           
           SEXP LECV
           ,
           SEXP pvalue,
           SEXP lower,
           SEXP give_log,
           SEXP PermutedStatistics
       )
       ;
extern 
       SEXP R_MaximumTest
       (
           /* R LECV Input */
           
           SEXP LECV
           ,
           SEXP alternative,
           SEXP pvalue,
           SEXP lower,
           SEXP give_log,
           SEXP PermutedStatistics,
           SEXP maxpts,
           SEXP releps,
           SEXP abseps
       )
       ;
extern 
       SEXP R_MaximallySelectedTest
       (
           SEXP LECV,
           SEXP ordered,
           SEXP teststat,
           SEXP minbucket,
           SEXP lower,
           SEXP give_log
       )
       ;
extern 
       SEXP R_ExpectationInfluence
       (
           /* R y Input */
           
               SEXP y,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       ) 
       ;
extern 
       SEXP R_CovarianceInfluence
       (
           /* R y Input */
           
               SEXP y,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           ,
           SEXP varonly
       )
       ;
extern 
       SEXP R_ExpectationX
       (
           /* R x Input */
           
               SEXP x,
           
           SEXP P,
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_CovarianceX
       (
           /* R x Input */
           
               SEXP x,
           
           SEXP P,
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           ,
           SEXP varonly
       )
       ;
extern 
       SEXP R_Sums
       (
           /* R N Input */
           
               SEXP N,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_KronSums
       (
           /* R x Input */
           
               SEXP x,
           
           SEXP P,
           /* R y Input */
           
               SEXP y,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           ,
           SEXP symmetric
       ) 
       ;
extern 
       SEXP R_KronSums_Permutation
       (
           /* R x Input */
           
               SEXP x,
           
           SEXP P,
           /* R y Input */
           
               SEXP y,
           
           /* R subset Input */
           
               SEXP subset
           ,
           SEXP subsety
       )
       ;
extern 
       SEXP R_colSums
       (
           /* R x Input */
           
               SEXP x,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_OneTableSums
       (
           /* R x Input */
           
               SEXP x,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_TwoTableSums
       (
           /* R x Input */
           
               SEXP x,
           
           /* R y Input */
           
               SEXP y,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_ThreeTableSums
       (
           /* R x Input */
           
               SEXP x,
           
           /* R y Input */
           
               SEXP y,
           
           /* R block Input */
           
               SEXP block
           ,
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           
       )
       ;
extern 
       SEXP R_order_subset_wrt_block
       (
           /* R y Input */
           
               SEXP y,
           
           /* R weights Input */
           
               SEXP weights
           ,
           /* R subset Input */
           
               SEXP subset
           ,
           /* R block Input */
           
               SEXP block
           
       )
       ;
extern 
       SEXP R_kronecker
       (
           SEXP A,
           SEXP B
       )
       ;

