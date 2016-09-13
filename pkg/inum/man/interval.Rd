\name{interval}
\alias{interval}
\alias{interval.numeric}
\title{
Cut Numeric Vectors in Intervals.
}
\description{
\code{interval} divides \code{x} into intervals and,
unlike \code{cut}, represents these as a numeric vector.
}
\usage{
interval(x, ...)
\method{interval}{numeric}(x, breaks = 50, ...)
}
\arguments{
  \item{x}{
A numeric vector.
}
  \item{breaks}{
  Either a numeric vector of two or more unique cut points or a
          single number (greater than or equal to 2) giving the number
          of intervals into which \code{x} is to be cut by \code{cut}.
}
  \item{\dots}{
Additional arguments, currently ignored.
}
}
\details{
This is just a wrapper around \code{cut} where the resulting intervals 
are stored as numeric values for simplified computation.
}
\value{
An object of class \code{interval}. A value of \code{0} encodes \code{NA}.
}
\seealso{
\code{\link{cut}}
}
\examples{
interval(0:100/100, breaks = 5)
}
\keyword{data}
