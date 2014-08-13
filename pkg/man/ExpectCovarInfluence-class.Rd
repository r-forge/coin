\name{ExpectCovarInfluence-class}
\docType{class}
\alias{ExpectCovarInfluence-class}
\title{Class \code{"ExpectCovarInfluence"}}
\description{
  Objects of class \code{"ExpectCovarInfluence"} represent the expectation and
  covariance of the influence function.
}
% NOTE: the markup in the following section is necessary for correct rendering
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \preformatted{
    new("ExpectCovarInfluence", \dots).}
}
\section{Slots}{
  \describe{
    \item{\code{sumweights}:}{
      Object of class \code{"numeric"}.  The sum of the case weights.
    }
    \item{\code{expectation}:}{
      Object of class \code{"numeric"}.  The expectation of the influence
      function.
    }
    \item{\code{covariance}:}{
      Object of class \code{"matrix"}.  The covariance matrix of the influence
      function.
    }
    \item{\code{dimension}:}{
      Object of class \code{"integer"}.  The dimension of the influence function.
    }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{ExpectCovar}"}, directly.
}
\section{Methods}{
  No methods defined with class \code{"ExpectCovarInfluence"} in the signature.
}
\section{Note}{
  Normally, this class is used internally only.
}
\keyword{classes}
