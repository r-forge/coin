\name{ScaleTests}
\alias{ansari_test.formula}
\alias{ansari_test.IndependenceProblem}
\alias{ansari_test}
\alias{fligner_test.formula}
\alias{fligner_test.IndependenceProblem}
\alias{fligner_test}
\title{ Independent Two- and K-Sample Scale Tests }
\description{
    Testing the equality of the distributions of a numeric response in
    two or more independent groups against scale alternatives.
}
\usage{

\method{ansari_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{ansari_test}{IndependenceProblem}(x, 
    alternative = c("two.sided", "less", "greater"),
    distribution = c("asympt", "approx", "exact"),
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...)

\method{fligner_test}{formula}(formula, data, subset = NULL, 
    weights = NULL, \dots)
\method{fligner_test}{IndependenceProblem}(x, 
    ties.method = c("mid-ranks", "average-scores"),
    distribution = c("asympt", "approx"), 
    ...)

}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where \code{y}
    is a numeric variable giving the data values and \code{x} a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
      integer valued weights for the observations.}
  \item{x}{an object of class \code{IndependenceProblem}.}
  \item{alternative}{the alternative hypothesis must be
    one of \code{"two.sided"} (default), \code{"greater"} or    
    \code{"less"}.  You can specify just the initial letter.}   
  \item{distribution}{the null distribution of the test statistic
    can be computed (\code{exact}) or can be approximated by its
    \code{asymptotic} distribution (\code{asympt}) 
    or via Monte-Carlo resampling (\code{approx}).}
  \item{ties.method}{two methods are available to adjust scores for ties.
    Either the score generating function is applied to \code{mid-ranks}
    or scores, based on random ranks, are averaged \code{average-scores}.}
  \item{conf.int}{a logical indicating whether a confidence interval
    for the difference in location should be computed.}
  \item{conf.level}{confidence level of the interval.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the distribution of \code{y} in
  the groups given by \code{x} is tested. In particular, the methods
  documented here are designed to detect scale alternatives. 

  For the adjustment of scores for tied values see Hajek, Sidak and Sen
  (1999), page 131ff.


}
\value{
  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue}, \code{statistic} and (optionally)
  \code{confint}.

}
\references{

        Bernd Streitberg & Joachim R\"ohmel (1986),
        Exact distributions for permutations and rank tests:
        An introduction to some recently published algorithms.
        \emph{Statistical Software Newsletter} \bold{12}(1), 10--17.

        Bernd Streitberg & Joachim R\"ohmel (1987),
        Exakte Verteilungen f\"ur Rang- und Randomisierungstests
        im allgemeinen $c$-Stichprobenfall.
        \emph{EDV in Medizin und Biologie} \bold{18}(1), 12--19.

        Mark A. van de Wiel (2001), The split-up algorithm: a fast
        symbolic method for computing p-values of rank statistics.
        \emph{Computational Statistics} \bold{16}, 519--538.

        David F. Bauer (1972), Constructing confidence sets using
        rank statistics. \emph{Journal of the American Statistical Association}
        \bold{67}, 687--690.

        Jaroslav Hajek, Zbynek Sidak & Pranab K. Sen (1999),
        \emph{Theory of Rank Tests}. San Diego, London: Academic Press.

}
\examples{

data(bloodp)

# Ansari-Bradley test

ansari_test(bp ~ group, data = bloodp)
ansari_test(bp ~ group, data = bloodp, distribution = "exact")
ansari_test(bp ~ group, data = bloodp, distribution = "approx", B = 9999)

# Lepage test, Hollander & Wolfe, page 169ff
platalet_counts <- data.frame(
    counts = c(120, 124, 215, 90, 67, 95, 190, 180, 135, 399, 
               12, 20, 112, 32, 60, 40),
    treatment = factor(c(rep("Prednisons", 10), rep("Control", 6))))

lt <- independence_test(counts ~ treatment, data = platalet_counts,
    ytrafo = function(data) trafo(data, numeric_trafo = function(x)       
        cbind(rank(x), ansari_trafo(x))),
    teststat = "quadtype", distribution = "approx", B = 9999)

lt

# where does the rejection come from?
lt2 <- independence_test(counts ~ treatment, data = platalet_counts,
    ytrafo = function(data) trafo(data, numeric_trafo = function(x) 
        matrix(c(rank(x), ansari_trafo(x)), ncol = 2,
               dimnames = list(1:length(x), c("Location", "Scale")))),
    teststat = "maxtype", distribution = "approx", B = 9999)

pvalue(lt2, adjusted = TRUE)

}
\keyword{htest}
