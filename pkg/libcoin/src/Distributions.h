
double C_chisq_pvalue(double stat, int df, int lower, int give_log);
double C_perm_pvalue(int greater, int B, int lower, int give_log);
double C_maxtype_pvalue(const double stat, const double *Covariance,
    const int n, int alternative, int lower, int give_log,
        int maxpts, double releps, double abseps, double tol);
void C_Permute(int *x, int n, int *ans);
void C_PermuteBlock(int *x, int *table, int Ntable, int *ans);
void C_doPermuteBlock(int *subset, int Nsubset, int *table, int Nlevels, 
                      int *Nsubset_tmp, int *perm);
void C_doPermute(int *subset, int Nsubset, int *Nsubset_tmp, int *perm);
void C_setup_subset(int N, int *N_ans);
void C_setup_subset_weights(int N, int *weights, int *sw_ans);
void C_setup_subset_weights_subset(int Nsubset, int *weights, int *subset, int *sw_ans);
void C_order_wrt_block(int *subset, int Nsubset, int *block, int *table, int Nlevels);
