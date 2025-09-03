\name{power.free1way.test}
\alias{rfree1way}
\alias{power.free1way.test}
\title{
Power Calculations for Distribution-free Wald Tests in Stratified One-Way Layouts
}
\description{
Compute power of tests via simulation and approximation or determine
parameters to obtain target power.
}
\usage{
rfree1way(n, prob = NULL, alloc_ratio = 1, blocks = ifelse(is.null(prob), 1, NCOL(prob)), 
          strata_ratio = 1, delta = 0, offset = 0, 
          link = c("logit", "probit", "cloglog", "loglog"))
power.free1way.test(n = NULL, prob = rep.int(1 / n, n), 
          alloc_ratio = 1, blocks = NCOL(prob),
          strata_ratio = 1, delta = NULL, mu = 0, 
          sig.level = .05, power = NULL,
          link = c("logit", "probit", "cloglog", "loglog"),
          alternative = c("two.sided", "less", "greater"), 
          nsim = 100, seed = NULL, tol = .Machine$double.eps^0.25) 
}
\arguments{
  \item{n}{number of observations in the control group of the first block.
}
  \item{prob}{a matrix defining the density of discrete controls for each
    block (in columns).
}
  \item{alloc_ratio}{allocation ratio, a factor defining the number of
    observations in each treatment group in the first block relative to
    \code{n}.
}
  \item{blocks}{number of blocks.
}
  \item{strata_ratio}{stratification ratio, a factor defining the number 
    of observations in each block relative to the first block.
}
  \item{delta}{true treatment effects comparing each treatment group to the
    control group.
}
  \item{offset}{a treatment-specific offset term, its length is recycled to
    the number of groups.
}
 \item{link}{a character defining a link function and thus the model and
    parameter interpretation. See \code{\link{free1way}}. 
}
\item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.
}
  \item{mu}{A vector indicating the true values of the treatment effect
    defining the null hypotheses.
}
  \item{sig.level}{significance level (Type I error probability).
}
  \item{power}{power of test (1 minus Type II error probability).
}
  \item{nsim}{number of simulations used to approximate the Hessian
    evaluated at the true treatment effects \code{delta}.
}
  \item{seed}{an object specifying if and how the random number generator
    should be initialized, see \code{\link[stats]{simulate}}.
}
  \item{tol}{numerical tolerance used in root finding, the default
    providing (at least) four significant digits.
}
}
\details{

The two functions use the same interface to sample from a specific
distribution-free semiparametric model (\code{rfree1way}) or to evaluate the
power of such a design (\code{power.free1way.test}). The latter function can
also be used to solve the power function for sample sizes or treatment
effects. The power function relies on an simulation-based approximation of the
Hessian. This is much faster than simulations for estimating the power
directly, but less accurate.

}
\value{
Object of class \code{power.htest}, a list of the arguments
(including the computed one) for \code{power.free1way.test}. A data
frame with simulated outcome values (\code{y}) for \code{rfree1way}.
}
\examples{

### sample from proportional odds model with 1:2 allocation
### based on odds ratio of 3
x <- rfree1way(n = 15, delta = log(3), alloc_ratio = 2)
### Wilcoxon-Mann-Whitney rank sum test
summary(free1way(y ~ groups, data = x), test = "Permutation")$p.value
wilcox.test(y ~ groups, data = x, exact = FALSE, correct = FALSE)$p.value

### approximate power for this test
power.free1way.test(n = 15, delta = log(3), alloc_ratio = 2)

}
\keyword{htest}
