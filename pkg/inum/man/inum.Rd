\name{inum}
\alias{inum}
\alias{inum.data.frame}
\title{
Coerse Variables in Data Frames to \code{enum} or \code{interval}
}
\description{
Represents elements of a data frame as \code{enum} or \code{interval}.
}
\usage{
inum(object, nmax = 20, ...)
\method{inum}{data.frame}(object, nmax = 20, ignore = NULL, 
     total = FALSE, weights = NULL, as.interval = "",
     complete.cases.only = FALSE, meanlevels = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A data frame.
}
  \item{nmax}{
Maximal number of categories for each of the numeric variables.
}
  \item{ignore}{
A character vector of variable names not to be discretised.
}
  \item{total}{
A logical. \code{TRUE} means that a condensed data frame of all 
variables is returned, \code{FALSE} a list of discretised variables.
}
  \item{weights}{
An optional vector of weights.
}
  \item{as.interval}{
A character vector of variable names to be converted to \code{\link{interval}}
instead of \code{\link{enum}}.
}
  \item{complete.cases.only}{
A logical. \code{TRUE} removes all rows with missing values.
}
  \item{meanlevels}{
A logical. \code{TRUE}, the level is the mean of the observations in the
corresponding bin. The default \code{FALSE} uses the largest observation in
the bin.
}
  \item{\dots}{
Additional arguments, currently ignored.
}
}
\details{
Each variable in \code{object} is converted to \code{\link{enum}} or \code{\link{interval}}.
}
\value{
An object of class \code{inum}, basically a list of
\code{\link{enum}} or \code{\link{interval}} objects. If
\code{total = TRUE}, an integer vector with a data frame
as \code{levels} attribute is returned. In this case, \code{0}
means \code{NA}.
}
\examples{
data("iris", package = "datasets")
iris[1,1] <- NA
inum(iris, nmax = 5)
inum(iris, nmax = 5, total = TRUE)
inum(iris, nmax = 5, total = TRUE, as.interval = "Sepal.Width",
     complete.cases.only = TRUE)
}
\keyword{data}
