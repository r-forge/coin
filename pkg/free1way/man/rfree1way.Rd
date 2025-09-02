\name{power.free1way.test}
\alias{rfree1way}
\alias{power.free1way.test}
\title{
Power Calculations for Distribution-free Wald Tests in Stratified One-Way Layouts
}
\description{
Compute power of test via simulation and approximation or determine
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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observations in the control group of the first block
}
  \item{prob}{A matrix defining the density of discrete controls for each
  block (in columns).
}
  \item{alloc_ratio}{Allocation ratio, a factor defining the number of
  observations in each treatment group in the first block relative to
  \code{n}
}
  \item{blocks}{Number of blocks
}
  \item{strata_ratio}{Stratification ratio, a factor defining the number 
  of observations in each block relative to the first block
}
  \item{delta}{True treatment effects comparing each treatment group to the
  control group
}
  \item{offset}{A treatment-specific offset term, its length is recycled to
  the number of groups
}
 \item{link}{A character defining a link function and thus the model and
    parameter interpretation. See \code{\link{free1way}}. }
\item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.}
  \item{mu}{A number indicating the true value of the treatment effect
  defining the null hypothesis
}
  \item{sig.level}{significance level (Type I error probability)
}
  \item{power}{power of test (1 minus Type II error probability)
}
  \item{nsim}{
%%     ~~Describe \code{link} here~~
}
  \item{seed}{
%%     ~~Describe \code{link} here~~
}
  \item{tol}{numerical tolerance used in root finding, the default
          providing (at least) four significant digits.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Object of class \code{power.htest}, a list of the arguments
(including the computed one) for \code{power.free1way.test}. A data
frame with simulated outcome values (\code{y}) for \code{rfree1way}.
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
