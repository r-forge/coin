\name{Distribution}
\alias{asymptotic}
\alias{approximate}
\alias{exact}
\title{Specification of the Reference Distribution}
\description{
  Specification of the asymptotic, approximate (Monte Carlo) and exact reference
  distribution.
}
\usage{
asymptotic(maxpts = 25000, abseps = 0.001, releps = 0)
approximate(B = 10000, parallel = c("no", "multicore", "snow"),
            ncpus = 1, cl = NULL)
exact(algorithm = c("auto", "shift", "split-up"), fact = NULL)
}
\arguments{
  \item{maxpts}{
    an integer, the maximum number of function values (see
    \code{\link[mvtnorm:algorithms]{GenzBretz}}).  Defaults to \code{25000}.
  }
  \item{abseps}{
    a numeric, the absolute error tolerance (see
    \code{\link[mvtnorm:algorithms]{GenzBretz}}).  Defaults to \code{0.001}.
  }
  \item{releps}{
    a numeric, the relative error tolerance (see
    \code{\link[mvtnorm:algorithms]{GenzBretz}}).  Defaults to \code{0}.
  }
  \item{B}{
    an integer, the number of Monte Carlo replicates used for the computation of
    the approximate reference distribution.  Defaults to \code{10000}.
  }
  \item{parallel}{
    a character, specifying the type of parallel operation: either \code{"no"}
    (default), \code{"multicore"} or \code{"snow"}.
  }
  \item{ncpus}{
    an integer, the number of processes to be used in parallel operation.
    Defaults to \code{1}.
  }
  \item{cl}{
    an object inheriting from class \code{"cluster"}, specifying an optional
    \pkg{parallel} or \pkg{snow} cluster if \code{parallel = "snow"}.  Defaults
    to \code{NULL}.
  }
  \item{algorithm}{
    a character, specifying the algorithm to be used for the computation of the
    exact reference distribution: either \code{"auto"} (default), \code{"shift"}
    or \code{"split-up"}.
  }
  \item{fact}{
    an integer to multiply the response values with.  Defaults to \code{NULL}.
  }
}
\details{
  The \code{distribution} argument to, e.g., \code{\link{independence_test}} can
  be specified with additional arguments using the functions documented here.

  The approximate (Monte Carlo) reference distribution is obtained through
  serial operation by default (\code{parallel = "no"}).  The use of parallel
  operation is specified by setting \code{parallel} to either \code{"multicore"}
  (not available for MS Windows) or \code{"snow"}.  In the latter case, if
  \code{cl = NULL} (default) a cluster with \code{ncpus} processes is created on
  the local machine unless a default cluster has been registered (see
  \code{\link[parallel:makeCluster]{setDefaultCluster}}) in which case that gets
  used instead.  Alternatively, the use of a particular \pkg{parallel} or
  \pkg{snow} cluster can be specified using \code{cl}.  (See \sQuote{Examples}.)
  For further details on the use of parallel operation, we refer to the
  \pkg{parallel} package and its documentation.

  Exact algorithms are currently implemented for two-sample problems only.
}
\examples{
## Approximate (Monte Carlo) Cochran-Mantel-Haenszel test

## Serial operation
set.seed(123)
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000))

\dontrun{
## Multicore with 8 processes (not for MS Windows)
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000,
                                    parallel = "multicore", ncpus = 8))

## Automatic PSOCK cluster with 4 processes
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000,
                                    parallel = "snow", ncpus = 4))

## Registered FORK cluster with 12 processes (not for MS Windows)
fork12 <- parallel::makeCluster(12, "FORK") # set-up cluster
parallel::setDefaultCluster(fork12) # register default cluster
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000,
                                    parallel = "snow"))
parallel::stopCluster(fork12) # clean-up

## User-specified PSOCK cluster with 8 processes
psock8 <- parallel::makeCluster(8, "PSOCK") # set-up cluster
set.seed(123, kind = "L'Ecuyer-CMRG")
cmh_test(disease ~ smoking | gender, data = alzheimer,
         distribution = approximate(B = 100000,
                                    parallel = "snow", cl = psock8))
parallel::stopCluster(psock8) # clean-up}
}
\keyword{htest}
