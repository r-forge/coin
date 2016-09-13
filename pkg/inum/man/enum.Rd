\name{enum}
\alias{enum}
\title{
  Enumeration-type Representation of Vectors
}
\description{
  Elements of a vector are stored as a set of levels and 
  an integer representing the enumeration.
}
\usage{
  enum(x)
}
\arguments{
  \item{x}{
  A vector. Currently, methods for factors, logicals, integers, and
  numeric vectors are implemented.
}
}
\details{
The unique elements of \code{x} are stored as a \code{levels} attribute
to an integer representing the enumeration. \code{levels} and \code{nlevels}
methods are available. This is essentially the same as \code{factor} where 
the levels can be arbitrary vectors, not just characters.
}
\value{
An object of class \code{enum}. A value of \code{0} encodes \code{NA}.
}
\seealso{
\code{\link{factor}}
}
\examples{

(ex <- enum(x <- gl(2, 2)))
all.equal(levels(ex)[ex], x)

(ex <- enum(x <- rep(c(TRUE, FALSE), 2)))
all.equal(levels(ex)[ex], x)

(ex <- enum(x <- rep(1:5, 2)))
all.equal(levels(ex)[ex], x)

(ex <- enum(x <- rep(1:5 + .5, 2)))
all.equal(levels(ex)[ex], x)

(ex <- enum(x <- c(NA, rep(1:5 + .5, 2))))
all.equal(c(NA, levels(ex))[unclass(ex) + 1L], x)

}
\keyword{data}
