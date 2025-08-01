\name{NullDistribution}
\alias{asymptotic}
\alias{approximate}
\alias{exact}
\title{Specification of the Reference Distribution}
\description{
  Specification of the asymptotic, approximative (Monte Carlo) and exact
  reference distribution.
}
\usage{
asymptotic(maxpts = 25000, abseps = 0.001, releps = 0)
approximate(nresample = 10000L, parallel = c("no", "multicore", "snow"),
            ncpus = 1L, cl = NULL, B)
exact(algorithm = c("auto", "shift", "split-up"), fact = NULL)
}
\arguments{
  \item{maxpts}{
    an integer, the maximum number of function values.  Defaults to
    \code{25000}.
  }
  \item{abseps}{
    a numeric, the absolute error tolerance.  Defaults to \code{0.001}.
  }
  \item{releps}{
    a numeric, the relative error tolerance.  Defaults to \code{0}.
  }
  \item{nresample}{
    a positive integer, the number of Monte Carlo replicates used for the
    computation of the approximative reference distribution.  Defaults to
    \code{10000L}.
  }
  \item{parallel}{
    a character, the type of parallel operation: either \code{"no"} (default),
    \code{"multicore"} or \code{"snow"}.
  }
  \item{ncpus}{
    an integer, the number of processes to be used in parallel operation.
    Defaults to \code{1L}.
  }
  \item{cl}{
    an object inheriting from class \code{"cluster"}, specifying an optional
    \pkg{parallel} or \pkg{snow} cluster if \code{parallel = "snow"}.  Defaults
    to \code{NULL}.
  }
  \item{B}{
    deprecated, use \code{nresample} instead.
  }
  \item{algorithm}{
    a character, the algorithm used for the computation of the exact reference
    distribution: either \code{"auto"} (default), \code{"shift"} or
    \code{"split-up"}.
  }
  \item{fact}{
    an integer to multiply the response values with.  Defaults to \code{NULL}.
  }
}
\details{
  \code{asymptotic()}, \code{approximate()} and \code{exact()} can be supplied
  to the \code{distribution} argument of, e.g.,
  \code{\link{independence_test}()} to provide control of the specification of
  the asymptotic, approximative (Monte Carlo) and exact reference distribution,
  respectively.

  The asymptotic reference distribution is computed using a randomised
  quasi-Monte Carlo method (Genz and Bretz, 2009) and is applicable to arbitrary
  covariance structures with dimensions up to 1000.  See
  \code{\link[mvtnorm:algorithms]{GenzBretz}()} in package \pkg{mvtnorm} for
  details on \code{maxpts}, \code{abseps} and \code{releps}.

  The approximative (Monte Carlo) reference distribution is obtained by a
  conditional Monte Carlo procedure, i.e., by computing the test statistic for
  \code{nresample} random samples from all admissible permutations of the
  response \eqn{\bf{Y}}{Y} within each block (Hothorn \emph{et al.}, 2008).  By
  default, the distribution is computed using serial operation
  (\code{parallel = "no"}).  The use of parallel operation is specified by
  setting \code{parallel} to either \code{"multicore"} (not available for MS
  Windows) or \code{"snow"}.  In the latter case, if \code{cl = NULL} (default)
  a cluster with \code{ncpus} processes is created on the local machine unless a
  default cluster has been registered (see
  \code{\link[parallel:makeCluster]{setDefaultCluster}()} in package
  \pkg{parallel}) in which case that gets used instead.  Alternatively, the use
  of an optional \pkg{parallel} or \pkg{snow} cluster can be specified by
  \code{cl}.  See \sQuote{Examples} and package \pkg{parallel} for details on
  parallel operation.

  The exact reference distribution, currently available for univariate
  two-sample problems only, is computed using either the shift algorithm
  (Streitberg and \enc{Röhmel}{Roehmel}, 1984, 1986, 1987) or the split-up
  algorithm (van de Wiel, 2001).  The shift algorithm handles blocks pertaining
  to, e.g., pre- and post-stratification, but can only be used with positive
  integer-valued scores \eqn{h(\bf{Y})}{h(Y)}.  The split-up algorithm can be
  used with non-integer scores, but does not handle blocks.  By default, an
  automatic choice is made (\code{algorithm = "auto"}) but the shift and
  split-up algorithms can be selected by setting \code{algorithm} to
  \code{"shift"} or \code{"split-up"}, respectively.
}
\note{
  Starting with version 1.1-0, the default for \code{algorithm} is
  \code{"auto"}, having identical behaviour to \code{"shift"} in previous
  versions.  In earlier versions of the package, \code{algorithm = "shift"}
  silently switched to the split-up algorithm if non-integer scores were
  detected, whereas the current version exits with a warning.

  In versions prior to 1.3-0, the number of Monte Carlo replicates in
  \code{approximate()} was specified using the now deprecated \code{B} argument.
  \strong{This will be made defunct and removed in a future release.}  It has
  been replaced by the \code{nresample} argument (for consistency with the
  \pkg{libcoin}, \pkg{party} and \pkg{partykit} packages).
}
\references{
  Genz, A. and Bretz, F.  (2009).  \emph{Computation of Multivariate Normal and
    t Probabilities}.  Heidelberg: Springer-Verlag.

  Hothorn, T., Hornik, K., van de Wiel, M. A. and Zeileis, A.  (2008).
  Implementing a class of permutation tests: The coin package.  \emph{Journal of
    Statistical Software} \bold{28}(8), 1--23.  \doi{10.18637/jss.v028.i08}

  Streitberg, B. and \enc{Röhmel}{Roehmel}, J.  (1984).  Exact nonparametrics
  in APL.  \emph{APL Quote Quad} \bold{14}(4), 313--325.
  \doi{10.1145/384283.801115}

  Streitberg, B. and \enc{Röhmel}{Roehmel}, J.  (1986).  Exact distributions
  for permutations and rank tests: an introduction to some recently published
  algorithms.  \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

  Streitberg, B. and \enc{Röhmel}{Roehmel}, J.  (1987).  Exakte verteilungen
  \enc{für}{fuer} rang- und randomisierungstests im allgemeinen
  c-stichprobenfall.  \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

  van de Wiel, M. A.  (2001).  The split-up algorithm: a fast symbolic method
  for computing p-values of distribution-free statistics.  \emph{Computational
    Statistics} \bold{16}(4), 519--538.  \doi{10.1007/s180-001-8328-6}
}
\examples{
## Approximative (Monte Carlo) Cochran-Mantel-Haenszel test

## Serial operation
set.seed(123)
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(nresample = 100000))

\dontrun{
## Multicore with 8 processes (not for MS Windows)
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(nresample = 100000,
                                    parallel = "multicore", ncpus = 8))

## Automatic PSOCK cluster with 4 processes
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(nresample = 100000,
                                    parallel = "snow", ncpus = 4))

## Registered FORK cluster with 12 processes (not for MS Windows)
fork12 <- parallel::makeCluster(12, "FORK") # set-up cluster
parallel::setDefaultCluster(fork12) # register default cluster
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(nresample = 100000,
                                    parallel = "snow"))
parallel::stopCluster(fork12) # clean-up

## User-specified PSOCK cluster with 8 processes
psock8 <- parallel::makeCluster(8, "PSOCK") # set-up cluster
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(nresample = 100000,
                                    parallel = "snow", cl = psock8))
parallel::stopCluster(psock8) # clean-up}
}
\keyword{htest}
