\name{free1way}
\alias{free1way}
\alias{free1way.formula}
\alias{free1way.table}
\alias{free1way.numeric}
\alias{free1way.factor}
\alias{print.free1way}
\alias{coef.free1way}
\alias{vcov.free1way}
\alias{logLik.free1way}
\alias{summary.free1way}
\alias{confint.free1way}
\title{
Distribution-free Inference in a Stratified One-Way Layout
}
\description{
Estimation, tests and confidence intervals for parameters in the
distribution-free stratified K-sample one-way layout for binary, ordinal,
numeric and potentially independently right-censored outcomes, including
semiparametrically efficient score tests against Lehmann, odds ratio, or
hazard ratio alternatives along with corresponding confidence intervals.
}
\usage{
free1way(y, ...)
\S3method{free1way}{formula}(formula, data, weights, subset, na.action = na.pass, ...)
\S3method{free1way}{numeric}(y, groups, blocks = NULL, event = NULL, weights = NULL, nbins = 0, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(groups)), 
                 deparse1(substitute(blocks))), ...)
\S3method{free1way}{factor}(y, groups, blocks = NULL, event = NULL, weights = NULL, 
    varnames = c(deparse1(substitute(y)), 
                 deparse1(substitute(groups)), 
                 deparse1(substitute(blocks))), ...)
\S3method{free1way}{table}(y, link = c("logit", "probit", "cloglog", "loglog"), 
    mu = 0, B = 0, ...)
\S3method{summary}{free1way}(object, test, alternative = c("two.sided", "less", "greater"), 
    tol = .Machine$double.eps, ...)
\S3method{coef}{free1way}(object, what = c("shift", "PI", "AUC", "OVL"), ...)
\S3method{confint}{free1way}(object, parm,
    level = .95, test = c("Permutation", "Wald", "LRT", "Rao"), 
    what = c("shift", "PI", "AUC", "OVL"), ...)
\S3method{vcov}{free1way}(object, ...)
\S3method{logLik}{free1way}(object, ...)
}
\arguments{
  \item{y}{a binary factor, an ordered factor, a numeric vector containing
    the response values or a \code{\link{table}} containing the response in the first
    dimension, the groups in the second dimension, and, optionally, strata and
    event indicators as third and fourth dimensions.}
  \item{nbins}{an optional integer defining the number of intervals to
    divide the range of a numeric response \code{y} into. The default is to
    cut the observations at breaks given by the uniquely observed values (\code{nbins = 0}). In
    the presence of right-censoring, uniquely observed event times define
    the breaks.}
  \item{groups}{a grouping factor.}
  \item{blocks}{a stratification factor, optional.}
  \item{event}{a logical vector representing events (\code{TRUE}) and
    independently right-censored observations (\code{FALSE}), optional.}
  \item{formula}{a formula of the form \code{y ~ x + strata(z)} where \code{y} gives 
    the sample values (binary, ordered or numeric) and \code{x} the corresponding groups.
    In stratified designs, a \code{strata(z)} term specifies the strata.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.  By default the variables are taken from
    \code{environment(formula)}.}
  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Should be \code{NULL} or a numeric vector.  If non-NULL,
    the weighted log-likelihood is maximised.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{varnames}{a character vector giving the names of the response,
    grouping, and stratifying varibales.}
  \item{link}{a character defining a link function and thus the model and
    parameter interpretation. See \sQuote{Details}.}
  \item{mu}{a number specifying an optional parameter used to form the
    null hypothesis.  See \sQuote{Details}.}	
  \item{B}{an integer specifying the number of replicates used in the
    permutation test. The default (\code{B = 0}) performs an asymptotic
    permutation test, \code{B > 0} uses \code{B} Monte Carlo replications to
    approximate the permutation distribution.}
  \item{object}{an object of class \code{free1way} as returned by
    \code{free1way}.}
  \item{test}{a character vector defining the global test procedure for all
    parameters: \code{"Permutation"} performs a conditional permutation score test
    under the randomisation model, \code{"Wald"} performs a Wald test,
    \code{"LRT"} a likelihood ratio, and \code{"Rao"} a Rao score test under the
    population model.}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{what}{a character defining a monotone transformation of the
    shift parameters: probabilistic indices (\code{"PI"}), area under the
    curves (\code{"AUC"}), or overlap coefficients (\code{"OVL"}). The
    default is to return parameters on the original shift scale.}
  \item{tol}{a positive numeric tolerance.}
  \item{parm}{a specification of which parameters are to be given
          confidence intervals, either a vector of numbers or a vector
          of names.  If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{\dots}{additional arguments to \code{\link{optim}}, such as a list of
    \code{control} parameters for the \code{method = "L-BFGS-B"}.}
}
\details{

The distribution function \eqn{F_1} of the response in the control group
(defined by the first factor level of the grouping factor \code{x}) is
compared to the distribution functions in the remaining groups \eqn{F_k} for
\eqn{k = 2, \dots, K}.  No assumptions are made regarding the form of \eqn{F_1},
all inference procedures are thus distribution-free.  However, a
semiparametric model assumes that each distribution function \eqn{F_k} only
depends on \eqn{F_1} and a scalar parameter \eqn{\beta_k}.  The link
function \eqn{Q} defines the model
\deqn{F_k(y) = Q^{-1}(Q(F_1(y)) - \beta_k)}
such that positive values of the parameter \eqn{\beta_k} correspond to
stochasically larger response values in group \eqn{k} when compared to the
control group.

The argument \code{link} is used to define the link function and thus
specific models: Log-odds ratio alternatives are based on the logit link
(\code{link = "logit"})
\deqn{\log\left(\frac{F_k(y)}{1 - F_k(y)}\right) = \log\left(\frac{F_1(y)}{1 - F_1(y)}\right) - \beta_k,}
log-hazard ratio alternatives on the complementary log-log link (\code{link =
"cloglog"})
\deqn{1 - F_k(y) = (1 - F_1(y))^{\exp(- \beta_k)},} Lehmann alternatives on
the log-log link (\code{link = "loglog"})
\deqn{F_k(y) = F_1(y)^{\exp(- \beta_k)},} and a shift alternative on a
latent normal scale on the probit link (\code{link = "probit"})
\deqn{\Phi^{-1}(F_1(Y_k)) \sim N(\beta_k, 1).}
If strata are present, the distribution functions may be stratum-specific
but the shift parameters are assumed to be constant across strata, such that
the above model holds for each stratum.

The control distribution function \eqn{F_1} is treated as a nuisance
parameter and the shift parameters \eqn{\beta_2, \dots, \beta_K} are
estimated by empirical maximum-likelihood estimation, the maximised
log-empirical likelihood is available via \code{logLik}. Parameter estimates and the
inverse observed Fisher information can be inspected via \code{coef} and
\code{vcov}.

Several global test procedures for the null hypothesis that all distribution
functions are identical (and thus \eqn{\beta_2 = \dots = \beta_K = 0}) 
can be specified by the \code{test} argument. Parameter-specific confidence
intervals are obtained via the inversion of a specific test procedure.

Parameter interpretation might be easier on transformed scales, such as
odds-ratios or hazard ratios. In addition, probabilistic indices (\code{what
= "PI"}, being equivalent to the area under the curve \code{"AUC"})
or overlap coefficients (\code{what = "OVL"}) are best obtained by inverting
permutation score, likelihood ratio, or Rao score test.

Assuming one of the semiparametric models, the parameter estimates are
semiparametrically efficient and the corresponding score tests and
confidence intervals are locally most powerful (Chapter 15.5 in van der
Vaart, 1998).

}
\value{
An object of class \code{free1way} with corresponding \code{logLik},
\code{coef}, \code{vcov}, \code{summary}, and \code{confint} methods.
}
\references{

  Aart W. Van der Vaart (1998). \emph{Asymptotic Statistics}. 
  Cambridge, UK: Cambridge University Press.

}
\examples{

### Friedman test
friedman.test(extra ~ group | ID, data = sleep)
ft <- free1way(extra ~ group | ID, data = sleep)
print(ft)
summary(ft)

}
\keyword{htest}
