\name{SurvTest}
\alias{surv_test}
\alias{logrank_test}
\alias{logrank_test.formula}
\alias{logrank_test.IndependenceProblem}
\title{ Independent Two- and K-Sample Tests for Censored Data }

\encoding{UTF-8}

\description{
  Testing the equality of survival distributions in two or more
  independent groups.
}
\usage{
\method{logrank_test}{formula}(formula, data, subset = NULL,
    weights = NULL, \dots)
\method{logrank_test}{IndependenceProblem}(object,
    ties.method = c("mid-ranks", "Hothorn-Lausen", "average-scores"),
    type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
             "Prentice-Marek", "Andersen-Borgan-Gill-Keiding",
             "Fleming-Harrington", "Self"),
    rho = NULL, gamma = NULL, \dots)
}
\arguments{
  \item{formula}{a formula of the form \code{Surv(time, event) ~ x | block} where
    \code{time} is a positive numeric variable denoting the survival time and
    \code{event} is a logical being \code{TRUE} when the event of interest
    was observed and \code{FALSE} in case of censoring. \code{x} is a factor
    with two or more levels giving the corresponding groups. \code{block} is an
    optional factor for stratification.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{weights}{an optional formula of the form \code{~ w} defining
    integer valued weights for the observations.}
  \item{object}{an object of class \code{IndependenceProblem}.}
  \item{ties.method}{a character specifying the way ties are handled in
    the definition of the logrank scores, see \sQuote{Details}.}
  \item{type}{a character specifying the type of test, see
    \sQuote{Details}.}
  \item{rho}{a numeric, the \eqn{\rho} constant when \code{type} is
    \code{"Tarone-Ware"}, \code{"Fleming-Harrington"} or
    \code{"Self"}.  Defaults to 0.5 in the first case and 0 in the other
    two cases.}
  \item{gamma}{a numeric, the \eqn{\gamma} constant when \code{type} is
    \code{"Fleming-Harrington"} or \code{"Self"}.  Defaults to 0.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{

  The null hypothesis of the equality of the survival distributions in
  the groups induced by \code{x} is tested.  The test implemented here
  is based on the weighted logrank test, but reformulated as a linear
  rank test (see \enc{Léton}{Leton} & Zuluaga, 2001, for the connection
  between the two formulations).  The family of weighted logrank tests
  encompasses a large collection of tests commonly used in the analysis
  of survival data.

  Let \eqn{(t_i, \delta_i)}, \eqn{i = 1, 2, \ldots, n}, represent a
  censored random sample of size \eqn{n}, where \eqn{t_i} is the
  observed survival time and \eqn{\delta_i} is the status indicator
  (\eqn{\delta_i} is 0 for censored observations and 1 otherwise).  To
  allow for ties in the data, let \eqn{t_{(1)} < t_{(2)} < \cdots <
    t_{(m)}} represent the \eqn{m}, \eqn{m \le n}, ordered distinct
  event times.  At time \eqn{t_{(k)}}, \eqn{k = 1, 2, \ldots, m}, the
  number of events and the number of subjects at risk are given by
  \eqn{d_k = \sum_{i = 1}^n I\!\left(t_i = t_{(k)}\,|\, \delta_i =
    1\right)} and \eqn{n_k = n - r_k} respectively, where \eqn{r_k}
  depends on the ties handling method.

  Three different methods of handling ties are available using
  \code{ties.method}: mid-ranks (\code{"mid-ranks"}, default), the
  Hothorn-Lausen method (\code{"Hothorn-Lausen"}) and average-scores
  (\code{"average-scores"}).  The first and last method are discussed
  and contrasted by Callaert (2003), whereas the second method is
  defined in Hothorn & Lausen (2003).  The mid-ranks method leads to
  \deqn{
    r_k^\mathrm{MR} = \sum_{i = 1}^n I\!\left(t_i < t_{(k)}\right)
  }
  and the Hothorn-Lausen method uses
  \deqn{
    r_k^\mathrm{HL} = \sum_{i = 1}^n I\!\left(t_i \le t_{(k)}\right) - 1
  .}
  The scores assigned to censored and uncensored observations at the
  \eqn{k}th event time are given by
  \deqn{
    C_k = \sum_{j = 1}^k w_j \frac{d_j}{n_j}
    \quad \mbox{and} \quad
    c_k = C_k - w_k
  }
  respectively, where \eqn{w} is the logrank weight.  In the
  average-scores method, used by, e.g., StatXact, the \eqn{d_k} events
  observed at the \eqn{k}th event time are arbitrarily ordered by
  assigning them distinct values \eqn{t_{(k_l)}}, \eqn{l = 1, 2, \ldots,
    d_k}, infinitesimally to the left of \eqn{t_{(k)}}.  Then scores
  \eqn{C_{k_l}} and \eqn{c_{k_l}} are computed as indicated above,
  effectively assuming that no event times are tied.  The scores
  \eqn{C_k} and \eqn{c_k} are assigned the average of the scores
  \eqn{C_{k_l}} and \eqn{c_{k_l}} respectively.

  It then follows that the score for the \eqn{i}th subject is
  \deqn{
    a_i = \left\{
      \begin{array}{ll}
        C_{k'} & \mbox{if } \delta_i = 0 \\
        c_{k'} & \mbox{otherwise}
      \end{array}
    \right.
  }
  where \eqn{k' = \max \{k: t_i \ge t_{(k)}\}}.  Thus, in the two-sample
  situation where \eqn{U_i = \{0, 1\}} denotes the group membership, the
  linear test statistic \eqn{T = \sum_{i = 1}^n a_i U_i} can be used to
  test the two-sided hypothesis \eqn{H_0\!: \theta = 1}, where
  \eqn{\theta = \lambda_2 / \lambda_1} and \eqn{\lambda_g} is the hazard
  rate in the \eqn{g}th group.  In case \code{alternative = "less"},
  i.e., the survival is lower in group 1 compared to group 2, the null
  hypothesis \eqn{H_0\!: \theta \ge 1} is tested.  The null hypothesis
  \eqn{H_0\!: \theta \le 1} is tested when
  \code{alternative = "greater"}, i.e., the survival is higher in group
  1 compared to group 2.

  The \code{type} argument allows for a choice between some of the most
  well-known members of the family of weighted logrank tests, each
  corresponding to a particular weight function.  The standard,
  unweighted, logrank test (\code{"logrank"}, default) was suggested by
  Mantel (1966), Peto & Peto (1972) and Cox (1972) and has \eqn{w_k =
    1}.  The Gehan-Breslow test (\code{"Gehan-Breslow"}) first proposed
  by Gehan (1965), and later extended to \eqn{K} samples by Breslow
  (1970), is a generalization of the Wilcoxon rank-sum test, where
  \eqn{w_k = n_k}.  The Tarone-Ware class of tests
  (\code{"Tarone-Ware"}) discussed by Tarone & Ware (1977) has \eqn{w_k
    = n_k^\rho}, where \eqn{\rho} is a constant; \eqn{\rho = 0.5}
  (default) was suggested by Tarone & Ware (1977), but note that
  \eqn{\rho = 0} and \eqn{\rho = 1} leads to the the standard logrank
  test and Gehan-Breslow test respectively.  The Prentice test
  (\code{"Prentice"}) is another generalization of the Wilcoxon rank-sum
  test suggested by Prentice (1978), where
  \deqn{
    w_k = \prod_{j = 1}^k \frac{n_j}{n_j + d_j}.
  }
  The Prentice-Marek test (\code{"Prentice-Marek"}) is yet another
  generalization of the Wilcoxon rank-sum test discussed by Prentice &
  Marek (1979), with
  \deqn{
    w_k = \prod_{j = 1}^k \frac{n_j + 1 - d_j}{n_j + 1}.
  }
  The Andersen-Borgan-Gill-Keiding test
  (\code{"Andersen-Borgan-Gill-Keiding"}) proposed by Andersen, Borgan,
  Gill & Keiding (1982) is a modified version of the Prentice-Marek test
  using
  \deqn{
    w_k = \frac{n_k}{n_k + 1}
            \prod_{j = 0}^{k - 1} \frac{n_j + 1 - d_j}{n_j + 1}
  }
  where \eqn{n_0 \equiv n} and \eqn{d_0 \equiv 0}.  The
  Fleming-Harrington class of tests (\code{"Fleming-Harrington"}
  suggested by Fleming & Harrington (1991) uses \eqn{w_k =
    \hat{S}_k^\rho (1 - \hat{S}_k)^\gamma}, where \eqn{\rho} and
  \eqn{\gamma} are constants and
  \deqn{
    \hat{S}_k = \prod_{j = 0}^{k - 1} \frac{n_j - d_j}{n_j},
    \quad
    \hat{S}_0 \equiv 1
  }
  is the \emph{left-continuous} Kaplan-Meier estimator of the survival
  function; \eqn{\rho = 0} and \eqn{\gamma = 0} leads to the standard
  logrank test.  The Self class of tests (\code{"Self"}) proposed by
  Self (1991) has \eqn{w_k = v_k^\rho (1 - v_k)^\gamma}, where
  \deqn{
    v_k = \frac{1}{2} \frac{t_{(k - 1)} + t_{(k)}}{t_{(m)}},
    \quad
    t_{(0)} \equiv 0
  }
  is the standardized mid-point between the \eqn{(k - 1)}th and the
  \eqn{k}th event time.  (This is a slight generalization of Self's
  original proposal in order to allow for non-integer follow-up times.)
  Again, \eqn{\rho} and \eqn{\gamma} are constants and \eqn{\rho = 0}
  and \eqn{\gamma = 0} leads to the standard logrank test.

  When \code{x} is an ordered factor, the linear-by-linear association
  test for ordered alternatives is computed.  In this case the test
  statistic is a scalar and the alternative hypothesis can be specified
  using the \code{alternative} argument.  The default scores,
  \code{1:nlevels(x)}, can be altered using the \code{scores} argument
  (see \code{\link{independence_test}}).  This procedure was first
  suggested by Tarone (1975) for the standard, unweighted, logrank test
  and later extended to general weights by Tarone & Ware (1977).

}
\note{

  The test statistic differs from the one used by
  \code{\link[survival]{survdiff}} (in the \pkg{survival} package), since the
  conditional variance is not identical to the variance estimate used by the
  classical logrank test.

  Combining \code{\link{independence_test}} or \code{\link{symmetry_test}} with
  \code{\link{logrank_trafo}} offers more flexibility than \code{logrank_test}
  and allows for, among other things, maximum-type versatile test procedures
  (e.g., Lee, 1996; see \sQuote{Examples}) and user-supplied logrank weights
  (see \code{\link{GTSG}} for tests against Weibull-type or crossing-curve
  alternatives).

  Starting with version 1.1-0, \code{logrank_test} replaces \code{surv_test},
  \strong{which is now deprecated and will be removed in a future release.}
  Furthermore, \code{logrank_trafo} is now an increasing function for all
  choices of \code{ties.method}, implying that the test statistic has the same
  sign irrespective of the ties handling method.  Consequently, the sign of the
  test statistic will now be the opposite of what it was in earlier versions
  unless \code{ties.method = "average-scores"}.  (In versions of \pkg{coin}
  prior to 1.1-0, \code{logrank_trafo} was a decreasing function of the data
  when \code{ties.method} was other than \code{"average-scores"}.)

}
\value{
  An object inheriting from class \code{"\linkS4class{IndependenceTest}"}.
}
\references{

  E. \enc{Léton}{Leton} & P. Zuluaga (2001).
  Equivalence Between Scores and Weighted Tests for Survival Curves.
  \emph{Communications in Statistics -- Theory and Methods}
  \bold{30}(4), 591--608.

  Herman Callaert (2003).
  Comparing Statistical Software Packages: The Case of the Logrank Test
  in StatXact.
  \emph{The American Statistician} \bold{57}(3), 214--217.

  Torsten Hothorn & Berthold Lausen (2003).
  On the Exact Distribution of Maximally Selected Rank Statistics.
  \emph{Computational Statistics & Data Analysis} \bold{43}(2),
  121--137.

  Nathan Mantel (1966).
  Evaluation of Survival Data and Two New Rank Order Statistics
  Arising in Its Consideration.
  \emph{Cancer Chemotherapy Reports} \bold{50}(3), 163--170.

  Richard Peto & Julian Peto (1972).
  Asymptotic Efficient Rank Invariant Test Procedures (with
  Discussion).
  \emph{Journal of the Royal Statistical Society} A \bold{135}(2),
  185--207.

  D. R. Cox (1972).
  Regression Models and Life-Tables (with Discussion).
  \emph{Journal of the Royal Statistical Society} B \bold{34}(2),
  187--220.

  Edmund A. Gehan (1965).
  A Generalized Wilcoxon Test for Comparing Arbitrarily
  Single-Censored Samples.
  \emph{Biometrika} \bold{52}(1--2), 203--223.

  Norman Breslow (1970).
  A Generalized Kruskal-Wallis Test for Comparing \eqn{K} Samples
  Subject to Unequal Patterns of Censorship.
  \emph{Biometrika} \bold{57}(3), 579--594.

  Robert E. Tarone & James Ware (1977).
  On Distribution-Free Tests for Equality of Survival Distributions.
  \emph{Biometrika} \bold{64}(1), 156--160.

  R. L. Prentice (1978).
  Linear Rank Tests with Right Censored Data.
  \emph{Biometrika} \bold{65}(1), 167--179.

  R. L. Prentice & P. Marek (1979).
  A Qualitative Descrepancy between Censored Data Rank Tests.
  \emph{Biometrics} \bold{35}(4), 861--867.

  Per Kragh Andersen, \enc{Ørnulf}{Ornulf} Borgan, Richard Gill & Niels
  Keiding (1982).
  Linear Nonparametric Tests for Comparison of Counting Processes,
  with Applications to Censored Survival Data (with Discussion).
  \emph{International Statistical Review} \bold{50}(3), 219--258.

  Thomas R. Fleming & David P. Harrington (1991).
  \emph{Counting Processes and Survival Analysis}.
  New York: John Wiley & Sons.

  Steven G. Self (1991).
  An Adaptive Weighted Log-Rank Test with Application to Cancer
  Prevention and Screening Trials.
  \emph{Biometrics} \bold{47}(3), 975--986.

  Robert E. Tarone (1975).
  Tests for Trend in Life Table Analysis.
  \emph{Biometrika} \bold{62}(3), 679--682.

  Jae Won Lee (1996).
  Some Versatile Tests Based on the Simultaneous Use of Weighted
  Log-Rank Statistics.
  \emph{Biometrics} \bold{52}(2), 721--725.

}
\examples{

  ### asymptotic tests for carcinoma data
  logrank_test(Surv(time, event) ~ stadium, data = ocarcinoma)
  survdiff(Surv(time, event) ~ stadium, data = ocarcinoma)


  ### example data given in Callaert (2003)
  exdata <- data.frame(time = c(1, 1, 5, 6, 6, 6, 6, 2, 2, 2, 3, 4, 4, 5, 5),
                       group = factor(c(rep(0, 7), rep(1, 8))))
  ### p = 0.0523
  survdiff(Surv(time) ~ group, data = exdata)
  ### p = 0.0505
  logrank_test(Surv(time) ~ group, data = exdata, distribution = "exact")
  ### p = 0.0468
  logrank_test(Surv(time) ~ group, data = exdata, distribution = "exact",
               ties.method = "average-scores")


  ### lung cancer example from StatXact
  lungcancer <-
      data.frame(time = c(257, 476, 355, 1779, 355,
                          191, 563, 242, 285, 16, 16, 16, 257, 16),
                 event = c(0, 0, 1, 1, 0,
                           1, 1, 1, 1, 1, 1, 1, 1, 1),
                 group = factor(rep(1:2, c(5, 9)),
                                labels = c("newdrug", "control")))

  ### StatXact 9 manual, page 214
  with(lungcancer,
       logrank_trafo(Surv(time, event), ties.method = "average-scores"))

  ### StatXact 9 manual, page 215
  logrank_test(Surv(time, event) ~ group, data = lungcancer,
               distribution = "exact", ties.method = "average-scores")

  ### StatXact 9 manual, page 222
  logrank_test(Surv(time, event) ~ group, data = lungcancer,
               distribution = "exact", ties.method = "average-scores",
               type = "Prentice")


  ### versatile test (Lee, 1996) using a range of rho and gamma
  rho.gamma <- expand.grid(rho = seq(-2, 2, 1), gamma = seq(0, 2, 1))
  vt <- independence_test(Surv(time, event) ~ group, data = lungcancer,
                          distribution = approximate(B = 10000),
                          ytrafo = function(data)
                              trafo(data, surv_trafo = function(y)
                                  logrank_trafo(y,
                                                ties.method = "average-scores",
                                                type = "Fleming-Harrington",
                                                rho = rho.gamma["rho"],
                                                gamma = rho.gamma["gamma"])))
  pvalue(vt, method = "step-down")

}
\keyword{htest}
