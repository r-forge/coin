\name{IndependenceProblem-class}
\docType{class}
\alias{IndependenceProblem-class}
\alias{coerce,IndependenceTestProblem,IndependenceProblem-method}
\alias{coerce<-,IndependenceTestProblem,IndependenceProblem-method}
\alias{initialize,IndependenceProblem-method}

\title{Class "IndependenceProblem" }
\description{  The main class for data handling. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("IndependenceProblem", x, y, block, weights)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"data.frame"}.}
    \item{\code{y}:}{Object of class \code{"data.frame"}.}
    \item{\code{weights}:}{Object of class \code{"numeric"}, however, 
      integer values are expected (and checked).}
    \item{\code{block}:}{Object of class \code{"factor"}.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "IndependenceTestProblem", to = "IndependenceProblem")}: ... }
    \item{coerce<-}{\code{signature(from = "IndependenceTestProblem", to = "IndependenceProblem")}: ... }
    \item{initialize}{\code{signature(.Object = "IndependenceProblem")}: ... }
  }
}
\keyword{classes}
