#include <R.h>
#include <Rdefines.h>

int psmirnov_exact_test_one(double q, double r, double s) {
    return ((r - s) >= q);
}

int psmirnov_exact_test_two(double q, double r, double s) {
    return (fabs(r - s) >= q);
}

double psmirnov_exact_uniq_lower(double q, int m, int n, int two) {
    double md, nd, *u, w;
    int i, j;
    int (*test)(double, double, double);

    md = (double) m;
    nd = (double) n;
    if(two)
	test = psmirnov_exact_test_two;
    else
	test = psmirnov_exact_test_one;

    u = (double *) R_alloc(n + 1, sizeof(double));

    u[0] = 1.;
    for(j = 1; j <= n; j++) {
        if(test(q, 0., j / nd))
            u[j] = 0.;
        else
            u[j] = u[j - 1];
    }
    for(i = 1; i <= m; i++) {
        w = (double)(i) / ((double)(i + n));
        if(test(q, i / md, 0.))
            u[0] = 0.;
        else
            u[0] = w * u[0];
        for(j = 1; j <= n; j++) {
            if(test(q, i / md, j / nd))
                u[j] = 0.;
            else
                u[j] = w * u[j] + u[j - 1];
        }
    }
    return u[n];
}

double psmirnov_exact_uniq_upper(double q, int m, int n, int two) {
    double md, nd, *u, v, w;
    int i, j;
    int (*test)(double, double, double);

    md = (double) m;
    nd = (double) n;
    if(two)
	test = psmirnov_exact_test_two;
    else
	test = psmirnov_exact_test_one;

    u = (double *) R_alloc(n + 1, sizeof(double));

    u[0] = 0.;
    for(j = 1; j <= n; j++) {
        if(test(q, 0., j / nd))
            u[j] = 1.;
        else
            u[j] = u[j - 1];
    }
    for(i = 1; i <= m; i++) {
        if(test(q, i / md, 0.))
            u[0] = 1.;
        for(j = 1; j <= n; j++) {
            if(test(q, i / md, j / nd))
                u[j] = 1.;
            else {
                v = (double)(i) / (double)(i + j);
                w = (double)(j) / (double)(i + j); /* 1 - v */
                u[j] = v * u[j] + w * u[j - 1];
            }
        }
    }
    return u[n];
}

double psmirnov_exact_ties_lower(double q, int m, int n, int *z, int two) {
    double md, nd, *u, w;
    int i, j;
    int (*test)(double, double, double);

    md = (double) m;
    nd = (double) n;
    if(two)
	test = psmirnov_exact_test_two;
    else
	test = psmirnov_exact_test_one;

    u = (double *) R_alloc(n + 1, sizeof(double));

    u[0] = 1.;
    for(j = 1; j <= n; j++) {
        if(test(q, 0., j / nd) && z[j])
            u[j] = 0.;
        else
            u[j] = u[j - 1];
    }
    for(i = 1; i <= m; i++) {
        w = (double)(i) / ((double)(i + n));
        if(test(q, i / md, 0.) && z[i])
            u[0] = 0.;
        else
            u[0] = w * u[0];
        for(j = 1; j <= n; j++) {
            if(test(q, i / md, j / nd) && z[i + j])
                u[j] = 0.;
            else
                u[j] = w * u[j] + u[j - 1];
        }
    }
    return u[n];
}

double psmirnov_exact_ties_upper(double q, int m, int n, int *z, int two) {
    double md, nd, *u, v, w;
    int i, j;
    int (*test)(double, double, double);

    md = (double) m;
    nd = (double) n;
    if(two)
	test = psmirnov_exact_test_two;
    else
	test = psmirnov_exact_test_one;

    u = (double *) R_alloc(n + 1, sizeof(double));

    u[0] = 0.;
    for(j = 1; j <= n; j++) {
        if(test(q, 0., j / nd) && z[j])
            u[j] = 1.;
        else
            u[j] = u[j - 1];
    }
    for(i = 1; i <= m; i++) {
        if(test(q, i / md, 0.) && z[i])
            u[0] = 1.;
        for(j = 1; j <= n; j++) {
            if(test(q, i / md, j / nd) && z[i + j])
                u[j] = 1.;
            else {
                v = (double)(i) / (double)(i + j);
                w = (double)(j) / (double)(i + j); /* 1 - v */
                u[j] = v * u[j] + w * u[j - 1];
            }
        }
    }
    return u[n];
}

SEXP psmirnov_exact(SEXP sq, SEXP sm, SEXP sn, SEXP sz,
		    SEXP stwo, SEXP slower) {
    double md, nd, p, q;
    int m, n, *z, two, lower;

    q = asReal(sq);
    m = asInteger(sm);
    n = asInteger(sn);
    two = asInteger(stwo);
    lower = asInteger(slower);

    md = (double) m;
    nd = (double) n;
    /*
       q has 0.5/mn added to ensure that rounding error doesn't
       turn an equality into an inequality, eg abs(1/2-4/5)>3/10 

    */
    q = (0.5 + floor(q * md * nd - 1e-7)) / (md * nd);

    if(sz == R_NilValue) {
	if(lower)
	    p = psmirnov_exact_uniq_lower(q, m, n, two);
	else
	    p = psmirnov_exact_uniq_upper(q, m, n, two);
    } else {
	z = INTEGER(sz);
	if(lower)
	    p = psmirnov_exact_ties_lower(q, m, n, z, two);
	else
	    p = psmirnov_exact_ties_upper(q, m, n, z, two);
    }
    return ScalarReal(p);
}
