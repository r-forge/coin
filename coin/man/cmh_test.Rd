\name{ContingencyTests}
\alias{cmh_test}
\alias{lbl_test}
\alias{chisq_test}
\title{ Independence in general I x K x J contingency tables }
\description{
    Testing the independence of two, possibly ordered, factors.
}
\usage{
\method{cmh_test}{formula}(formula, data, subset = NULL, \dots)
\method{cmh_test}{table}(x, distribution = c("asympt", "approx"), ...)
\method{cmh_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx"), ...)
\method{chisq_test}{formula}(formula, data, subset = NULL, \dots)
\method{chisq_test}{table}(x, distribution = c("asympt", "approx"), ...)
\method{chisq_test}{IndependenceProblem}(x, 
    distribution = c("asympt", "approx"), ...)
\method{lbl_test}{formula}(formula, data, subset = NULL, \dots)
\method{lbl_test}{table}(x, distribution = c("asympt", "approx"), ...)
\method{lbl_test}{IndependenceProblem}(x, 
    xscores = NULL, yscores = NULL, 
    distribution = c("asympt", "approx"), ...)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    and \code{x} are factors (possibly ordered) and \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{x}{an object inheriting from class \code{IndependenceProblem} or an
    object of class \code{table}.}
  \item{distribution}{the null distribution of the test statistic
    can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{xscores}{if \code{x} is ordered, a vector of length
      \code{nlevels(x)} with alternative scores.}
  \item{yscores}{if \code{y} is ordered, a vector of length
      \code{nlevels(y)} with alternative scores.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the independence of \code{y} and \code{x} is
  tested.

}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\examples{

data(jobsatisfaction)

### for females only
chisq_test(as.table(jobsatisfaction[,,"Female"]), 
    distribution = "approx", B = 9999)

### both Income and Job.Satisfaction unordered
cmh_test(jobsatisfaction)

### both Income and Job.Satisfaction ordered
lbl_test(jobsatisfaction)

### dose response data, StatXact 6 manual, 984 
tox <- as.table(matrix(c(100, 18, 50, 50, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1),
                       nrow = 4,
                       dimnames = list(dose = paste((1:4)*100, "mg"),
                                       tox = c("mild", "moderate", "severe", "death"))))

### 0.0708, page 993
lbl_test(tox)

### 0.0792, page 993
pvalue(lbl_test(tox, distribution = "approx", B = 9999))

### 0.0828, page 997
lbl_test(tox, scores = list(tox = c(1, 3, 9, 27)))

## 0.0780, page 997
pvalue(lbl_test(tox, scores = list(tox = c(1, 3, 9, 27)), distribution = "approx", B = 9999))

}
\keyword{htest}
