\name{coin}
\alias{coin}
\title{General Information on the coin Package}
\non_function{}
\description{

The \code{coin} package implements a general framework for conditional
inference procedures, commonly known as \emph{permutation tests},
theoretically derived by Strasser & Weber (1999). The conditional
expectation and covariance for a broad class of multivariate linear
statistics as well as the corresponding multivariate limiting distribution
was derived by Strasser & Weber (1999) which can be utilized to construct
tests for independence between two sets of variables. 

Beside a general implementation of the abstract framework the package offers
a rather huge set of convenience functions implementing well known classical
as well as less well known classical and non-classical test procedures
fitting in this framework. Representatives are linear rank statistics for
the two- and K-sample location and scale problem against ordered and
unordered alternatives including post-hoc tests for arbitrary contrasts,
tests of independence for contingency tables, two- and K-sample tests for
censored data, tests for independence of two continuous variables as well as 
tests for marginal homogenity and symmetry. Most of the classical procedures
given in famous text books like Hollander & Wolfe (1999) or Agresti (2002)
can be implemented as part of the general framework with much effort.
Approximations of the exact null-distribution via the limiting distribution and 
conditional Monte-Carlo procedures are available for every test while the 
exact null-distribution is currently available for two-sample problems only. 

}
\references{

    Helmut Strasser & Christian Weber (1999),
    On the asymptotic theory of permutation statistics. 
    \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.  

    Myles Hollander & Douglas A. Wolfe (1999),  
    \emph{Nonparametric Statistical Methods, 2nd Edition}.
    New York: John Wiley & Sons.

    Alan Agresti (2002), \emph{Categorical Data Analysis}. Hoboken, New
    Jersey: John Wiley & Sons.

}
\authors{

    Torsten Hothorn <Torsten.Hothorn@R-project.org>,
    Kurt Hornik <Kurt.Hornik@R-project.org>, 
    Mark van de Wiel <m.a.v.d.wiel@TUE.nl> and 
    Achim Zeileis <Achim.Zeileis@R-project.org>

}
\keyword{misc}
