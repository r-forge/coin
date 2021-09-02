
Pimp stats::ks.test

 - exact conditional p-values for two-sample Smirnov test in presence of
   ties

 - psmirnov (high-level R implementation of Schröer-Trenkler algorithm)
   replaces stats::C_pSmirnov2x

 - modifications to ks.test:
   - compute exact conditional p-values when n.x * n.y < 10000 (with or
     without ties) or when exact = TRUE
   - fire a warning when asymptotic p-values are reported in the presence of
     ties (asymptotic distribution assumes absolute continuous distribution)
   - use correct names of tests:
     one-sample => Kolmogorov-Smirnov test
     two-sample => Smirnov test
   - Report "Exact" or "Approximate" as part of the name of the test (this
     can be debated)
   - add formula interface and print data.name / alternative in terms
     of variable/group level names

 - modifications to ks.test.Rd:
   - use correct names (Kolmogorov-Smirnov vs Smirnov) and cite
     Berger & Zou overview paper
   - Ties handling: Differentiate between one-sample and two-sample
     situation. Explain exact conditional p-values for the two-sample
     case and cite Schröer and Trenkler (1995).
   - Add example with ties.

 - new confband generic / method for ks.test
   - implement nonparametric quantile-quantile plot and corresponding 
     confidence bands following Switzer (1976)
   - confband.Rd has two examples

Suggested modifications to stats:

 - remove C_pSmirnov2x from stats/src/ks.c (only covers the unconditional
   case and is only marginally faster compared to psmirnov())
   OR
   retain C_pSmirnov2x and call it in psmirnov when exact & is.null(obs)
   (for testing purposes, for example)
 - add psmirnov() to stats/R/ks.R (the original APL code is contained
   as comment and it would be good to keep this piece of code because the
   original source in the diploma thesis is almost impossible to obtain).
   Remove ::: (1x)
 - replace/add ks.test, ks.test.default and ks.test formula 
   in stats/R/ks.R after removing ::: (1x)
 - replace manual page stats/man/ks.test.Rd
 - update stats NAMESPACE

Optional additions:

 - add qsmirnov to stats/R/ks.R (quantile function)
 - add confband generic (complementing stats::confint generic)
 - add confband.ks.test and corresponding plot method (quantile-quantile
    plot with confidence bands)
