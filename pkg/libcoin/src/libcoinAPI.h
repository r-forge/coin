
#include <R_ext/Rdynload.h>
#include <libcoin.h>

// external API
void libcoin_R_LinstatExpCox(const SEXP data, const SEXP inputs,
                             const SEXP y, const SEXP weights) {

    static void(*fun)(SEXP, SEXP, SEXP, SEXP) = NULL;

    if (fun == NULL)
        fun = (void(*)(SEXP, SEXP, SEXP, SEXP)) R_GetCCallable("libcoin", "R_LinstatExpCov");

    fun(data, inputs, y, weights);
}
                                                                                                                                                                                                                                                                
