
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

 - modifications to ks.test.Rd:
   - use correct names (Kolmogorov-Smirnov vs Smirnov) and cite
     Berger & Zou overview paper
   - Ties handling: Differentiate between one-sample and two-sample
     situation. Explain exact conditional p-values for the two-sample
     case and cite Schröer and Trenkler (1995).
   - Add example with ties.

Suggested modifications to stats:

 - remove C_pSmirnov2x from stats/src/ks.c (only covers the unconditional
   case and is only marginally faster compared to psmirnov())
 - add psmirnov() to stats/R/ks.R
 - replace ks.test() in stats/R/ks.R
 - replace manual page stats/man/ks.test.Rd

