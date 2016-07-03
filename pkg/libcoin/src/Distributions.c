
#include "libcoin_internal.h"
#include "Tables.h"
#include "Sums.h"
#include "Utils.h"
#include <mvtnormAPI.h>

double C_chisq_pvalue(double stat, int df, int lower, int give_log)
{
    return(pchisq(stat, (double) df, lower, give_log));
}

/**
    Conditional asymptotic P-value of a maxabs-type test statistic\n
    Basically the functionality from package `mvtnorm' \n
    *\param stat test statitstic
    *\param Covariance covariance matrix
    *\param n nrow(Covariance)
    *\param maxpts number of Monte-Carlo steps
    *\param releps relative error
    *\param abseps absolute error
    *\param tol tolerance
*/

double C_maxabsstat_pvalue(const double stat, const double *Covariance, 
    const int n, int maxpts, double releps, double abseps, double tol) {

    int nu = 0, inform, i, j, sub, nonzero, *infin, *index, rnd = 0;
    double ans, myerror, *lower, *upper, *delta, *corr, *sd;

    /* univariate problem */
    if (n == 1) 
        return(2*pnorm(fabs(stat)*-1.0, 0.0, 1.0, 1, 0)); /* return P-value */
    
    if (n == 2)  
         corr = Calloc(1, double);
    else 
         corr = Calloc(n + ((n - 2) * (n - 1))/2, double);
    
    sd = Calloc(n, double);
    lower = Calloc(n, double);
    upper = Calloc(n, double);
    infin = Calloc(n, int);
    delta = Calloc(n, double);
    index = Calloc(n, int);

    /* determine elements with non-zero variance */ 

    nonzero = 0;
    for (i = 0; i < n; i++) {
        if (Covariance[S(i, i, n)] > tol) {
            index[nonzero] = i;
            nonzero++;
        }
    }

    /* mvtdst assumes the unique elements of the triangular 
       covariance matrix to be passes as argument CORREL 
    */

    for (int iz = 0; iz < nonzero; iz++) {

        /* handle elements with non-zero variance only */
        i = index[iz];

        /* standard deviations */
        sd[i] = sqrt(Covariance[S(i, i, n)]);
                
        /* always look at the two-sided problem */           
        lower[iz] = fabs(stat) * -1.0;
        upper[iz] = fabs(stat);
        infin[iz] = 2;
        delta[iz] = 0.0;
        
        /* set up vector of correlations, i.e., the upper 
           triangular part of the covariance matrix) */
        for (int jz = 0; jz < iz; jz++) {
            j = index[jz]; 
            sub = (int) (jz + 1) + (double) ((iz - 1) * iz) / 2 - 1;
            if (sd[i] == 0.0 || sd[j] == 0.0) 
                corr[sub] = 0.0; 
            else 
                corr[sub] = Covariance[S(i, j, n)] / (sd[i] * sd[j]);
        }
    }
        
    /* call mvtnorm's mvtdst C function defined in mvtnorm/include/mvtnormAPI.h */
    mvtnorm_C_mvtdst(&nonzero, &nu, lower, upper, infin, corr, delta, 
                     &maxpts, &abseps, &releps, &myerror, &ans, &inform, &rnd);

    /* inform == 0 means: everything is OK */
    switch (inform) {
        case 0: break;
        case 1: warning("cmvnorm: completion with ERROR > EPS"); break;
        case 2: warning("cmvnorm: N > 1000 or N < 1"); 
                ans = 0.0; 
                break;
        case 3: warning("cmvnorm: correlation matrix not positive semi-definite"); 
                ans = 0.0; 
                break;
        default: warning("cmvnorm: unknown problem in MVTDST");
                ans = 0.0;
    }
    Free(corr); Free(sd); Free(lower); Free(upper); 
    Free(infin); Free(delta); 
    return(1 - ans);  /* return P-value */
}


void C_Permute(int *x, int n, int *ans) 
{
    int k = n, j;
    
    for (int i = 0; i < k; i++) {
        j = n * unif_rand();
        ans[i] = x[j];
        x[j] = x[--n];
    }
}

void C_PermuteBlock(int *x, int *table, int Ntable, int *ans)
{
    int *px, *pans;
    
    px = x;
    pans = ans;
    
    for (int j = 0; j < Ntable; j++) {
        if (table[j] > 0) {
            C_Permute(px, table[j], pans);
            px += table[j];
            pans += table[j];
        }
    }
}

void C_doPermuteBlock(int *subset, int Nsubset, int *table, int Nlevels, 
                      int *Nsubset_tmp, int *perm) 
{
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_PermuteBlock(Nsubset_tmp, table, Nlevels, perm);
}

void C_doPermute(int *subset, int Nsubset, int *Nsubset_tmp, int *perm) 
{
    Memcpy(Nsubset_tmp, subset, Nsubset);
    C_Permute(Nsubset_tmp, Nsubset, perm);
}

void C_setup_subset(int N, int *N_ans)
{
    for (int i = 0; i < N; i++) N_ans[i] = i;
}

void C_setup_subset_weights(int N, int *weights, int *sw_ans)
{
    int itmp = 0;
    /* subset = rep(1:length(weights), weights) */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < weights[i]; j++)
            sw_ans[itmp++] = i;
    }
}

void C_setup_subset_weights_subset(int Nsubset, int *weights, int *subset, int *sw_ans)
{
    /* subset = rep(subset, weights[subset]) */
    int itmp = 0;
    for (int i = 0; i < Nsubset; i++) {
        for (int j = 0; j < weights[subset[i]]; j++)
            sw_ans[itmp++] = subset[i];
    }
}

void C_order_wrt_block(int *subset, int Nsubset, int *block, int *table, int Nlevels)
{
    int *cumtable, *subset_tmp;
    
    cumtable = Calloc(Nlevels, int);
    for (int k = 0; k < Nlevels; k++) cumtable[k] = 0;

    subset_tmp = Calloc(Nsubset, int);
    Memcpy(subset_tmp, subset, Nsubset);

    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++)
        cumtable[k] = cumtable[k - 1] + table[k - 1];
    
    for (int i = 0; i < Nsubset; i++)
        subset[cumtable[block[subset_tmp[i]]]++] = subset_tmp[i];
    Free(cumtable); Free(subset_tmp);
} 
