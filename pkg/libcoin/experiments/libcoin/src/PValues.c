
#include "libcoin.h"

double C_chisq_pvalue(double stat, int df, int give_log)
{
    return(pchisq(stat, df, 0, give_log));
}

        