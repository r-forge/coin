\name{pvalue-methods}
\docType{methods}
\alias{pvalue-methods}
\alias{pvalue,NullDistribution-method}
\alias{pvalue,IndependenceTest-method}
\alias{pvalue,ScalarIndependenceTest-method}
\alias{pvalue,MaxTypeIndependenceTest-method}
\alias{pvalue,QuadTypeIndependenceTest-method}
\title{ Extract p-values from test objects }
\description{

  Extracts the p-value from objects representing null distributions of 
  independence tests.

}
\section{Methods}{
\describe{

\item{object = "NullDistribution"}{ calls the \code{pvalue} slot from
    \code{object}.}

\item{object = "IndependenceTest"}{ calls the \code{pvalue} slot from
    \code{object@distribution}.}

\item{object = "ScalarIndependenceTest"}{ calls the \code{pvalue} slot from
    \code{object@distribution}.}

\item{object = "MaxTypeIndependenceTest"}{ calls the \code{pvalue} slot from
    \code{object@distribution}.}

\item{object = "QuadTypeIndependenceTest"}{ calls the \code{pvalue} slot from
    \code{object@distribution}.}

}}
\keyword{methods}
\keyword{htest}
