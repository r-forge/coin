
#include "coin_common.h"

/* Helpers.c */
/*
extern SEXP R_kronecker
(
    SEXP A,
    SEXP B
);
*/
extern SEXP R_maxstattrafo
(
    SEXP x,
    SEXP cutpoints
);

extern SEXP R_outersum
(
    SEXP A,
    SEXP B
);


/* StreitbergRoehmel.c */
extern SEXP R_cpermdist2
(
    SEXP score_a,
    SEXP score_b,
    SEXP m_a,
    SEXP m_b,
    SEXP retProb
);

extern SEXP R_cpermdist1
(
    SEXP scores
);


/* vandeWiel.c */
extern SEXP R_split_up_2sample
(
    SEXP scores,
    SEXP m,
    SEXP obs,
    SEXP tol
);
