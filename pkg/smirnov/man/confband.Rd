\name{confband}
\title{Confidence Bands}
\alias{confband}
\alias{confband.ks.test}
\description{
  Computes confidence bands for functions defined by a model.
  There is a method for objects inheriting from class
  \code{"\link{ks.test}"}.
}
\usage{
confband(object, level = 0.95, \dots)
}
\arguments{
  \item{object}{a fitted model object.}
  \item{level}{the confidence level required.}
  \item{\dots}{additional argument(s) for methods.}
}
\value{
  A matrix with columns giving the estimated function (\code{"Estimate"}) as
  well as lower (\code{"lwr"}) and upper (\code{"upr"}) confidence limits 
  for each function argument (\code{"x"}). Confidence bands control 
  the coverage probability (\code{"level"}) for all intervals 
  simultaneously.
}
\details{
  \code{confband} is a generic function.  

  For objects of class \code{"ks.test"} a confidence band for the treatment
  function in two-sample parallel designs is computed following Switzer
  (1976). This function is also known as quantile-quantile plot.

}
\references{

  Paul Switzer (1976). Confidence Procedures for Two-sample Problems.
  Biometrika \bold{63}(1), 13--25. \doi{10.1093/biomet/63.1.13}.

}
\seealso{
  \code{\link{ks.test}}
}
\examples{
# data from Switzer (1976)
angle <- c(-31, -30, -25, -25, -23, -23, -22, -20, -20, -18,
           -18, -18, -16, -15, -15, -14, -13, -11, -10, - 9,
           - 8, - 7, - 7, - 7, - 6, - 6, - 4, - 4, - 3, - 2,
           - 2, - 1,   1,   1,   4,   5,  11,  12,  16,  34,
           -31, -20, -18, -16, -16, -16, -15, -14, -14, -14,
           -14, -13, -13, -11, -11, -10, - 9, - 9, - 8, - 7,
           - 7, - 6, - 6,  -5, - 5, - 5, - 4, - 2, - 2, - 2,
             0,   0,   1,   1,   2,   4,   5,   5,   6,  17)
sex <- gl(2, 40, labels = c("Female", "Male"))

(sk <- ks.test(angle ~ sex))
cb <- confband(sk, level = .945)
plot(cb) # => no indication of deviation from the null because 
         # diagonal is covered by the confidence band

set.seed(29)
# artificial example
grp <- rep(gl(2, 1, labels = c("Group A", "Group B")), c(40, 60))
# shift model
x1 <- rnorm(length(grp), mean = c(1, 2)[grp], sd = 1)
# shift/scale model
x2 <- rnorm(length(grp), mean = c(1, 2)[grp], sd = c(1, .5)[grp])
layout(matrix(1:4, nrow = 2))
# treatment function is parallel to diagonal, `significant' 
# deviation from the null of both distributions being equal
# because diagonal not covered by confidence band
plot(confband(ks.test(x1 ~ grp)), main = "Shift model")
# true treatment function in blue
lines(-3:3, -3:3 + 1, col = "blue")
boxplot(x1 ~ grp)
# shift model not appropriate because treatment function
# is nonlinear, also `significant' deviation from null
plot(confband(ks.test(x2 ~ grp)), main = "Shift/scale model")
# true treatment function
lines(-3:3, qnorm(pnorm(-3:3, mean = 1, sd = 1), mean = 2, sd = .5), 
      col = "blue")
boxplot(x2 ~ grp)

}
\keyword{models}

