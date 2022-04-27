/**
    Exact Distribution of Two-Sample Permutation Tests
    Streitberg & Roehmel Algorithm

*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

/**
    Permutation density of the permutation two-sample test statistic

        sum(scores[1:m])

    for independent two sample problem. 

    REFERENCES

    Bernd Streitberg & Joachim R\"ohmel (1986),
    Exact distributions for permutations and rank tests:
    An introduction to some recently published algorithms.
    Statistical Software Newsletter 12(1), 10-17.

    Bernd Streitberg & Joachim R\"ohmel (1987),
    Exakte Verteilungen f\"ur Rang- und Randomisierungstests
    im allgemeinen $c$-Stichprobenfall.
    EDV in Medizin und Biologie 18(1), 12-19 (in german).

    *\param score (integer): increasingly ordered positive
                             score vector (typically ranks)
    *\param m (integer):     number of observations in 
                             the group defining the test statistic
*/

SEXP C_StreitbergRoehmel(SEXP score, SEXP m) {

    SEXP ret;                   /* density (return value) */
    SEXP iscore, im;
    int n, im_a, im_b = 0;      /* number of observations */
    int i, j, k, sum_a = 0, sum_b = 0, s_a = 0, s_b = 0, isb, idx, min_b, sum_b1;
    double msum = 0.0;          /* little helpers */
    int *iscore_b;              /* pointers */
    double *dH, *dret;

    iscore = PROTECT(coerceVector(score, INTSXP));
    n = LENGTH(iscore);
    /* check if sum(1:n) > INT_MAX */
    if (n > (-1 + sqrt(1 + 4 * 2 * (double)(INT_MAX - 1))) / 2)
        error("length of score vector too large");

    im = PROTECT(coerceVector(m, INTSXP));
    im_a = INTEGER(im)[0];

    iscore_b = INTEGER(iscore);

    /* compute the total sum of the scores and check if they are >= 0 */

    for (i = 0; i < n; i++) {
        if (iscore_b[i] < 0)
            error("score for observation number %d is negative", i);
        if (i > 0) {
            if (iscore_b[i] < iscore_b[i - 1])
                error("score is not ordered increasingly");
        }
        if (sum_b + iscore_b[i] > INT_MAX)
            error("sum of scores is larger than INT_MAX");
        sum_b += iscore_b[i];
    }
    sum_a = n * (n + 1) / 2;
    for (i = n - im_a; i < n; i++) 
        im_b += iscore_b[i];

    /*
      optimization according to Streitberg & Roehmel
    */

    sum_a = imin2(sum_a, im_a);
    sum_b = imin2(sum_b, im_b);
    sum_b1 = sum_b + 1;

    if ((double) (sum_a + 1) * (double) (sum_b + 1) > INT_MAX)
        error("dimension too large");

    /*
        initialize H
    */

    dH = Calloc((sum_a + 1) * (sum_b + 1), double);

    for (i = 0; i <= sum_a; i++) {
        isb = i * sum_b1;
        for (j = 0; j <= sum_b; j++) 
            dH[isb + j] = 0.0;
    }

    /*
        start the Shift-Algorithm with H[0,0] = 1
    */

    dH[0] = 1.0;

    for (k = 0; k < n; k++) {
        s_a += 1;
        s_b += iscore_b[k];

        /*
            compute H up to row im_a and column im_b
        */

        min_b = imin2(im_b,s_b);

        for (i = imin2(im_a, s_a); i >= 1; i--) {
            isb = i * sum_b1;
            idx = (i - 1) * sum_b1;
            for (j = min_b; j >= iscore_b[k]; j--)
                dH[isb + j] += dH[idx + (j - iscore_b[k])];
        }
    }

    PROTECT(ret = allocVector(REALSXP, sum_b));
    dret = REAL(ret);

    /*
        get the values for sample size im_a (in row m) and sum it up
    */

    isb = im_a * (sum_b + 1);
    for (j = 0; j < sum_b; j++) {
        if (!R_FINITE(dH[isb + j + 1]))
            error("overflow error; cannot compute exact distribution");
        dret[j] = dH[isb + j + 1];
        msum += dret[j];
    }
    if (!R_FINITE(msum) || msum == 0.0)
        error("overflow error; cannot compute exact distribution");

    /*
        compute probabilities and return the density 
    */

    for (j = 0; j < sum_b; j++)
        dret[j] = dret[j] / msum;

    Free(dH);

    UNPROTECT(3);

    return(ret);
}
