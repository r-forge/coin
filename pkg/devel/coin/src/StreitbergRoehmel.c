/**
    Exact Distribution of Two-Sample Permutation Tests
    Streitberg and Röhmel Shift Algorithm

    *\file StreitbergRoehmel.c
    *\author $Author$
    *\date $Date$
*/

#include "coin_common.h"
#include <Rmath.h> /* for imin2 */

/**
    The density of the permutation distribution for independent two-sample
    problems

    REFERENCES

    Streitberg, B. and Röhmel, J.  (1986).  Exact distributions for permutations
    and rank tests: An introduction to some recently published algorithms.
    Statistical Software Newsletter 12(1), 10-17.

    Streitberg, B. and Röhmel, J.  (1987).  Exakte Verteilungen für Rang- und
    Randomisierungstests im allgemeinen c-Stichprobenfall.  EDV in Medizin und
    Biologie 18(1), 12-19.

    *\param score_b score vector (typically ranks)
    *\param m_a integer indicating the sum of m_a elements of score_a
    *\param m_b integer indicating the sum of m_b elements of score_b
*/

SEXP R_cpermdist2(SEXP score_b, SEXP m_a) {
    /* Compute the joint permutation distribution of the sum of the first 'm_a'
       elements of 'score_a' and 'score_b'.  In this case the exact conditional
       distribution in the independent two-sample problem is computed. */

    /* number of observations */
    int n, im_a, im_b = 0;
    /* matrix of permutations and vector of probs */
    SEXP x;
    /* little helpers */
    int sum_a = 0, sum_b = 0, sum_bp1, s_a = 0, s_b = 0, isb;
    double msum = 0.0;
    /* pointers to R structures */
    int *iscore_b;
    double *dH, *dx;

    /* some basic checks, should be improved */
    if (!isVector(score_b))
        error("score_b is not a vector");

    n = LENGTH(score_b);

    iscore_b = INTEGER(score_b);

    im_a = INTEGER(m_a)[0];
    for (int i = n - im_a; i < n; i++)
        im_b += iscore_b[i];

    /* compute the total sum of the scores and check if they are >= 0 */
    sum_a = n; /* remember: score_a = (1,...,1) */
    for (int i = 0; i < n; i++) {
        if (iscore_b[i] < 0)
            error("score_b for observation number %d is negative", i);
        sum_b += iscore_b[i];
    }

    /* optimization according to Streitberg and Röhmel */
    sum_a = imin2(sum_a, im_a);
    sum_b = imin2(sum_b, im_b);

    /* initialize H */
    sum_bp1 = sum_b + 1;
    dH = R_Calloc((sum_a + 1) * sum_bp1, double);
    for (int i = 0; i <= sum_a; i++) {
        isb = i * sum_bp1;
        for (int j = 0; j <= sum_b; j++)
            dH[isb + j] = 0.0;
    }

    /* start the shift algorithm with H[0,0] = 1 */
    dH[0] = 1.0;
    for (int k = 0; k < n; k++) {
        s_a += 1; /* remember: score_a = (1,...,1) */
        s_b += iscore_b[k];
        /* compute H up to row im_a and column im_b
           Note: sum_a = min(sum_a, m) and sum_b = min(sum_b, c) */
        for (int i = imin2(im_a, s_a); i >= 1; i--) {
            isb = i * sum_bp1;
            for (int j = imin2(im_b, s_b); j >= iscore_b[k]; j--)
                dH[isb + j] +=
                    dH[(i - 1) * sum_bp1 + (j - iscore_b[k])];
        }
    }

    PROTECT(x = allocVector(REALSXP, sum_b));
    dx = REAL(x);
    /* get the values for sample size im_a (in row m) and sum it up */
    isb = im_a * (sum_b + 1);
    for (int j = 0; j < sum_b; j++) {
        if (!R_FINITE(dH[isb + j + 1]))
            error("overflow error; cannot compute exact distribution");
        dx[j] = dH[isb + j + 1];
        msum += dx[j];
    }
    if (!R_FINITE(msum) || msum == 0.0)
        error("overflow error; cannot compute exact distribution");
    /* compute probabilities and return the density x to R
       the support is min(score_b):sum(score_b) */
    for (int j = 0; j < sum_b; j++)
        dx[j] = dx[j] / msum;

    R_Free(dH);

    UNPROTECT(1);
    return(x);
}

/**
    The density of the permutation distribution for one sample problems.

    REFERENCES

    Streitberg, B. and Röhmel, J.  (1986).  Exact distributions for permutations
    and rank tests: An introduction to some recently published algorithms.
    Statistical Software Newsletter 12(1), 10-17.

    Streitberg, B. and Röhmel, J.  (1987).  Exakte Verteilungen für Rang- und
    Randomisierungstests im allgemeinen c-Stichprobenfall.  EDV in Medizin und
    Biologie 18(1), 12-19.

    *\param scores score vector (such as rank(abs(y)) for wilcoxsign_test)
*/

SEXP R_cpermdist1(SEXP scores) {
    /* compute the permutation distribution of the sum of the absolute values of
       the positive elements of 'scores' */

    /* number of observations */
    int n;
    /* vector giving the density of statistics 0:sum(scores) */
    SEXP H;
    /* little helpers */
    int sum_a = 0, s_a = 0;
    double msum = 0.0;
    /* pointers to R structures */
    int *iscores;
    double *dH;

    n = LENGTH(scores);
    iscores = INTEGER(scores);

    for (int i = 0; i < n; i++)
        sum_a += iscores[i];

    /* initialize H */
    PROTECT(H = allocVector(REALSXP, sum_a + 1));
    dH = REAL(H);
    for (int i = 0; i <= sum_a; i++)
        dH[i] = 0.0;

    /* start the shift algorithm with H[0] = 1.0 */
    dH[0] = 1.0;
    for (int k = 0; k < n; k++) {
        s_a = s_a + iscores[k];
            for (int i = s_a; i >= iscores[k]; i--)
                dH[i] = dH[i] + dH[i - iscores[k]];
    }

    /* get the number of permutations */
    for (int i = 0; i <= sum_a; i++) {
        if (!R_FINITE(dH[i]))
            error("overflow error: cannot compute exact distribution");
        msum += dH[i];
    }
    if (!R_FINITE(msum) || msum == 0.0)
        error("overflow error: cannot compute exact distribution");

    /* compute probabilities and return the density H to R */
    for (int i = 0; i <= sum_a; i++)
        dH[i] = dH[i] / msum; /* 0 is a possible realization */

    UNPROTECT(1);
    return(H);
}
