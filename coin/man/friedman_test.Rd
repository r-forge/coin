\name{SymmetryTests}
\alias{friedman_test}
\alias{mh_test}
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
\method{mh_test}{formula}(formula, data, subset = NULL, \dots)
\method{mh_test}{table}(x, yscores = NULL, ...)
\method{mh_test}{SymmetryProblem}(x, 
    distribution = c("asympt", "approx"), ...) 
\method{wilcoxsign_test}{formula}(formula, data, subset = NULL, \dots)
\method{wilcoxsign_test}{SymmetryProblem}(x, 
    distribution = c("asympt", "approx"), ...) 
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric (\code{friedman_test,wilcoxsign_test}) variable 
    or factor (\code{mh_test}) giving the data values and 
    \code{x} a factor with two (\code{wilcoxsign_test}) 
    or more levels giving the corresponding groups. \code{block} is an
    optional factor (which is generated automatically when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{x}{an object inheriting from class \code{SymmetryProblem} or a table 
           with identical dimnames.}
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
### Wilcoxon-Nemenyi-McDonald-Thompson
if (require(multcomp)) {

    ### setup SymmetryProblem (need to rank here)
    r <- as.vector(apply(matrix(RoundingTimes[["times"]], 
                                ncol = 3, byrow = TRUE),
                         1, rank))
    sp <- new("SymmetryProblem", RoundingTimes["methods"],
              data.frame(r = r))

    ### all pair comparisons
    rtt <- symmetry_test(sp, teststat = "maxtype",
         xtrafo = function(data)
             trafo(data, factor_trafo = function(x)
                 model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey"))
             )
     )

    ### a global test, again
    print(pvalue(rtt))

    ### simultaneous P-values for all pair comparisons
    ### Wide Angle vs. Round Out differ (Hollander and Wolfe, 1999, page 296)
    print(pvalue(rtt, adjusted = TRUE))
}

### strength of cotton, StatXact 6 manual, page 467
sc <- data.frame(block = factor(c(rep(1, 5), rep(2, 5), rep(3, 5))),
                 potash = ordered(rep(c(144, 108, 72, 54, 36), 3)),
                 strength = c(7.46, 7.17, 7.76, 8.14, 7.63,
                              7.68, 7.57, 7.73, 8.15, 8.00,
                              7.21, 7.80, 7.74, 7.87, 7.93))

### Page test
friedman_test(strength ~ potash | block, data = sc)

### approximate null distribution via Monte-Carlo
pvalue(friedman_test(strength ~ potash | block, data = sc, 
                     distribution = "approx", B = 9999))

### Opinions on Pre- and Extramarital Sex, Agresti (2000), page 421
opinions <- c("always wrong", "almost always wrong", 
              "wrong only sometimes", "not wrong at all")

PreExSex <- as.table(matrix(c(144, 33, 84, 126, 
                                2,  4, 14,  29, 
                                0,  2,  6,  25, 
                                0,  0,  1,  5), nrow = 4, 
                            dimnames = list(PremaritalSex = opinions,
                                            ExtramaritalSex = opinions)))

### treating response as nominal
mh_test(PreExSex)

### and as ordinal
mh_test(PreExSex, scores = list(response = 1:length(opinions)))

}
\keyword{htest}
