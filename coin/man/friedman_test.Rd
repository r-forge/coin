\name{SymmetryTests}
\alias{friedman_test}
\alias{bowker_test}
\alias{wilcoxsign_test}
\title{ Symmetry Tests }
\description{
    Testing the symmetry of a response for repeated measurements in a 
    complete block design.
}
\usage{
\method{friedman_test}{formula}(formula, data, subset = NULL, \dots)
\method{friedman_test}{SymmetryProblem}(x, 
    distribution = c("asympt", "approx"), ...) 
\method{bowker_test}{formula}(formula, data, subset = NULL, \dots)
\method{bowker_test}{SymmetryProblem}(x, 
    distribution = c("asympt", "approx"), ...) 
\method{wilcoxsign_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcoxsign_test}{SymmetryProblem}(x, 
    distribution = c("asympt", "approx"), ...) 
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric (\code{friedman_test,wilcoxsign_test}) variable 
    or factor (\code{bowker_test}) giving the data values and 
    \code{x} a factor with two (\code{wilcoxsign_test}) 
    or more levels giving the corresponding groups. \code{block} is an
    optional factor (which is generated automatically when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{x}{an object inheriting from class \code{SymmetryProblem}.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of symmetry of \code{y} across \code{x} is 
  tested.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\examples{

RoundingTimes <- data.frame(times = c(5.40, 5.50, 5.55,
    5.85, 5.70, 5.75,
    5.20, 5.60, 5.50,
    5.55, 5.50, 5.40,
    5.90, 5.85, 5.70,
    5.45, 5.55, 5.60,
    5.40, 5.40, 5.35,
    5.45, 5.50, 5.35,
    5.25, 5.15, 5.00,
    5.85, 5.80, 5.70,
    5.25, 5.20, 5.10,
    5.65, 5.55, 5.45,
    5.60, 5.35, 5.45,
    5.05, 5.00, 4.95,
    5.50, 5.50, 5.40,
    5.45, 5.55, 5.50,
    5.55, 5.55, 5.35,
    5.45, 5.50, 5.55,
    5.50, 5.45, 5.25,
    5.65, 5.60, 5.40,
    5.70, 5.65, 5.55,
    6.30, 6.30, 6.25),
    methods = factor(rep(c("Round Out", "Narrow Angle", "Wide Angle"), 22)))

### classical global test
friedman_test(times ~ methods, data = RoundingTimes)

### parallel coordinates plot
matplot(t(matrix(RoundingTimes$time, ncol = 3, byrow = TRUE)), 
        type = "l", col = 1, lty = 1, axes = FALSE, ylab = "Time", 
        xlim = c(0.5, 3.5))
axis(1, at = 1:3, labels = levels(RoundingTimes$methods))
axis(2)

### where do the differences come from?
if (require(multcomp)) {

    ### setup SymmetryProblem
    sp <- new("SymmetryProblem", RoundingTimes["methods"],
              RoundingTimes["times"])

    ### asymptotic Tukey all pair comparisons (without rank transformation)
    rtt <- symmetry_test(sp, teststat = "maxtype",
        xtrafo = function(data) 
            trafo(data, factor_trafo = function(x)
                model.matrix(~ x) \%*\% t(contrMat(table(x), "Tukey"))
            )
    )

    ### a global test, again
    pvalue(rtt)

    ### simultaneous P-values for all pair comparisons
    ### indicate _no_ difference between Narrow Angle and Round Out
    pvalue(rtt, adjusted = TRUE)
}

}
\keyword{htest}
