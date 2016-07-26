\name{doTest}
\alias{doTest}
\title{
  Permutation Test
}
\description{
  Perform permutation test for a linear statistic
}
\usage{
doTest(object, type = c("maxstat", "quadform"), 
       alternative = c("two.sided", "less", "greater"), 
       lower = FALSE, log = FALSE, 
       minbucket = 10L, ordered = TRUE, pargs = GenzBretz())
}
\arguments{
  \item{object}{an object returned by \code{\link{LinStatExpCov}}.
}
  \item{type}{use maximum-type or quadratic form as test statistic.
}
  \item{alternative}{alternative for maximum-type statistics.
}
  \item{lower}{a logical indicating if a p-value (\code{lower} is \code{FALSE})
               or 1 - p-value (\code{lower} is \code{TRUE}) shall be returned.
}
  \item{log}{a logical, if \code{TRUE} probabilities are log-probabilities.
}
  \item{minbucket}{minimum weight in either of two groups for maximally selected
                   statistics.
}
  \item{ordered}{a logical, if \code{TRUE} maximally selected statistics 
                 assume that the cutpoints are ordered.
}
  \item{pargs}{arguments as in \code{\link[mvtnorm]{GenzBretz}}.
}
}
\details{
  Computes a test statistic, a corresponding p-value and, optionally, cutpoints for 
  maximally selected statistics.
}
\value{
  A list.
}
\keyword{htest}