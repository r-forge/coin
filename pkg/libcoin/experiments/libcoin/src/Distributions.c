
#include "libcoin.h"
#include "Tables.h"
#include "Sums.h"
#include "helpers.h"

double C_chisq_pvalue(double stat, int df, int give_log)
{
    return(pchisq(stat, df, 0, give_log));
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
            sw_ans[++itmp] = subset[i];
    }
}

void C_order_wrt_block(int *subset, int Nsubset, int *block, int *table, int Nlevels)
{
    int *cumtable, *subset_tmp;
    
    cumtable = Calloc(Nlevels, int);
    subset_tmp = Calloc(Nsubset, int);
    Memcpy(subset_tmp, subset, Nsubset);
    cumtable[0] = 0;
    /* table[0] are missings, ie block == 0 ! */
    for (int k = 1; k < Nlevels; k++) cumtable[k] = cumtable[k - 1] + table[k - 1];
    
    for (int i = 0; i < Nsubset; i++)
        subset[cumtable[block[subset[i]]]++] = subset_tmp[i];
    Free(cumtable); Free(subset_tmp);
} 
