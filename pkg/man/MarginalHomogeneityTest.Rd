\name{MarginalHomogeneityTest}
\alias{mh_test}
\alias{mh_test.formula}
\alias{mh_test.table}
\alias{mh_test.SymmetryProblem}

\title{ Marginal Homogeneity Test }

\description{
  Testing marginal homogeneity in a complete block design.
}
\usage{
\method{mh_test}{formula}(formula, data, subset = NULL, \dots)
\method{mh_test}{table}(object, \dots)
\method{mh_test}{SymmetryProblem}(object, \dots)
}
\arguments{
  \item{formula}{a formula of the form \code{y ~ x | block} where
    \code{y} and \code{x} are factors (possibly ordered) giving the
    response and corresponding measurement conditions respectively.
    \code{block} is an optional factor (which is generated automatically
    when omitted).}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{object}{an object inheriting from class \code{"SymmetryProblem"}
    or an object of class \code{"table"} with identical \code{dimnames}
    attributes.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of marginal homogeneity for matched pairs or
  matched sets is tested.

  This procedure is known as the McNemar test (McNemar, 1947) when both
  \code{x} and \code{y} are binary factors.  When \code{x} is a binary
  factor and \code{y} is a factor with an arbitrary number of levels,
  the procedure is known as the Stuart(-Maxwell) test (Stuart, 1955;
  Maxwell, 1970); see also Bhapkar (1966) for a more powerful Wald-type
  test.  (As noted by Ireland, Ku & Kullback (1969), Bhapkar's statistic
  is given by \eqn{Q / (1 - Q / n)}, where \eqn{Q} is Stuart's statistic
  and \eqn{n} is the sample size).  If both \code{x} and \code{y} are
  factors with an arbitrary number of levels, the procedure is known as
  the Madansky test of interchangeability (Madansky, 1963), which
  implies marginal homogeneity.

  When both \code{x} and \code{y} are ordered factors, or can be
  regarded as such, the linear-by-linear association test for ordered
  alternatives is computed.  In this case the test statistic is a scalar
  and the alternative hypothesis can be specified using the
  \code{alternative} argument.  The default scores, \code{1:nlevels(x)}
  and \code{1:nlevels(y)} respectively, can be altered using the
  \code{scores} argument (see \code{\link{symmetry_test}}).  This
  extension was given by Birch (1965) who also discussed the situation
  when either the response or the measurement conditions are ordered
  factors; see also White, Landis & Cooper (1982).

  An overview of the these procedures is given by Agresti (2002).

}
\note{

  This function is rather inefficient for data with a large number of
  observations.

}
\value{

  An object inheriting from class \code{IndependenceTest} with
  methods \code{show}, \code{pvalue} and \code{statistic}.

}
\references{

  Quinn McNemar (1947).
  Note on the sampling error of the difference between correlated
  proportions or percentages.
  \emph{Psychometrika} \bold{12}(2), 153--157.

  Alan Stuart (1955).
  A test for homogeneity of the marginal distributions in a two-way
  classification.
  \emph{Biometrika} \bold{42}(3/4), 412--416.

  A. E. Maxwell (1970).
  Comparing the classification of subjects by two independent judges.
  \emph{British Journal of Psychiatry} \bold{116}(535), 651--655.

  V. P. Bhapkar (1966).
  A note on the equivalence of two test criteria for hypotheses in
  categorical data.
  \emph{Journal of the American Statistical Association} \bold{61}(313),
  228--235.

  C. T. Ireland, H. H. Ku & S. Kullback (1969).
  Symmetry and marginal homogeneity of an \eqn{r \times r} contingency
  table.
  \emph{Journal of the American Statistical Association} \bold{64}(328),
  1323--1341.

  Albert Madansky (1963).
  Tests of homogeneity for correlated samples.
  \emph{Journal of the American Statistical Association} \bold{58}(301),
  97--119.

  M. W. Birch (1965).
  The detection of partial association, II: the general case.
  \emph{Journal of the Royal Statistical Society} B \bold{27}(1),
  111--124.

  Andrew A. White, J. Richard Landis & Murray M. Cooper (1982).
  A note on the equivalence of several marginal homogeneity test
  criteria for categorical data.
  \emph{International Statistical Review} \bold{50}(1), 27--34.

  Alan Agresti (2002).
  \emph{Categorical Data Analysis}, 2nd Edition.
  Hoboken, New Jersey: John Wiley & Sons.

}
\examples{
  ### Performance of prime minister
  ### Agresti (2002, p. 409)
  rating  <- c("Approve", "Disprove")
  performance <- matrix(c(794, 150,
                           86, 570),
                        nrow = 2, byrow = TRUE,
                        dimnames = list(First = rating,
                                        Second = rating))

  ### Asymptotic McNemar Test
  diag(performance) <- 0 # speed-up: only off-diagonal elements contribute
  mh_test(as.table(performance))

  ### Exact McNemar Test
  mh_test(as.table(performance), distribution = "exact")


  ### Opinions on Pre- and Extramarital Sex, Agresti (2002), page 421
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
\dontshow{\dontrun{
  ### example taken from
  ### http://www.john-uebersax.com/stat/mcnemar.htm#stuart
  rating <- c("low", "moderate", "high")
  x <- as.table(matrix(c(20, 10,  5,
                         3, 30, 15,
                         0,  5, 40),
                       ncol = 3, byrow = TRUE,
                       dimnames = list(Rater1 = rating, Rater2 = rating)))
  ### test statistic W_0 = 13.76
  mh_test(x)
}}

  ### Vote intention example from Madansky (1963, p. 107--108)
  intention <- c("Republican","Democratic","Uncertain")
  vote <- as.table(array(c(120, 1,  8, 2,   2,  1, 2, 1,  7,
                             6, 2,  1, 1, 103,  5, 1, 4,  8,
                            20, 3, 31, 1,   6, 30, 2, 1, 81),
                         dim = c(3, 3, 3),
                         dimnames = list(July = intention,
                                         August = intention,
                                         June = intention)))
  ### test statistic Q = 70.77
  mh_test(vote)


  ### Cross-over study
  ### http://www.nesug.org/proceedings/nesug00/st/st9005.pdf
  relief <- c("none", "moderate", "complete")
  dysmenorrhea <- as.table(array(c(6, 2, 1,  4, 3, 0,  5, 2, 2,
                                   3, 1, 0, 13, 3, 0, 10, 1, 0,
                                   1, 2, 1,  8, 1, 1, 14, 2, 0),
                                 dim = c(3, 3, 3),
                                 dimnames =  list("Placebo" = relief,
                                                  "High dose" = relief,
                                                  "Low dose" = relief)))

  ### Ordered response, Q = 53.76
  mh_test(dysmenorrhea, scores = list(response = 1:3))

  ### Both response and measurement conditions ordered, Q = 47.29
  mh_test(dysmenorrhea, scores = list(response = 1:3, conditions = 3:1))

}
\keyword{htest}
