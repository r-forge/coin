\name{ppplot}
\alias{ppplot}
\title{
Probability-probability Plots
}
\description{
\code{ppplot} produces a probabiliy-probability (PP) plot of two numerical variables, essentially
the Receiver Operating Characteristic (ROC) curve. 
If \code{conf.level} is given, an estimate and corresponding 
confidence band for the ROC curve under a 
distribution-free semiparametric model is plotted.
}
\usage{
ppplot(x, y, plot.it = TRUE, xlab = deparse1(substitute(x)), 
       ylab = deparse1(substitute(y)), ..., conf.level = NULL, 
       conf.args = list(link = "logit", type = "Wald", col = NA, border = NULL))
}
\arguments{
  \item{x}{the first sample for \code{ppplot}, a numerical variable.
}
  \item{y}{the second sample for \code{ppplot}, a numerical variable.
}
  \item{plot.it}{logical. Should the result be plotted?
}
  \item{xlab, ylab}{the \code{xlab} and \code{ylab} refer to the y
          and x axes respectively.
}
  \item{\dots}{graphical parameters.
}
  \item{conf.level}{confidence level of the band. The default, \code{NULL}, does not
          lead to the computation of a confidence band.
}
  \item{conf.args}{list of arguments defining confidence band computation and
         visualisation: \code{link} defines the link function of a
         distribition-free semiparametric model, \code{type} specifies the
         statistical concept the confidence band is derived from; see
         \code{\link{free1way}} for other options. The remaining
         elements govern how the band is plotted.
}
}
\details{

The plot produces a ROC curve, however avoiding rearrangements pushing the curve above the
diagonal: The dots represent pairs \eqn{(\hat{F}_{x}(y), \hat{F}_{y}(y))}.

If the data generating process follows a model where the two distribution
functions, after appropriate transformation, are horizontally shifted
versions of each other, the ROC curve is a simple function of this shift and
confidence bands can be obtained from a confidence interval for this shift
parameter, see \code{free1way} for the model and
\bibcitet{free1way::SewakHothorn2023} for the connection to ROC curves.

Substantial deviations of the empirical (points) from the theoretical
(lines) ROC curve indicates lack of fit of the semiparametric model.

}
\value{
A list with components \code{x} and \code{y} representing the empirical ROC
curve.
}
\references{
  \bibshow{*}
}
\examples{

### well-fitting logistic model
nd <- data.frame(groups = gl(2, 50, labels = paste0("G", 1:2)))
nd$y <- rlogis(nrow(nd), location = c(0, 2)[nd$groups])
with(with(nd, split(y, groups)),
     ppplot(G1, G2, conf.level = .95,
            conf.args = list(link = "logit", type = "Wald", col = 2)))
## with appropriate Wilcoxon test and log-odds ratio
coef(ft <- free1way(y ~ groups, data = nd))
## area under the (model-based!) ROC curve (AUC) and confidence interval
coef(ft, what = "AUC")
confint(ft, test = "Rao", what = "AUC")
## the model-based ROC curve
prb <- 1:99 / 100
points(prb,  plogis(qlogis(prb) - coef(ft)), pch = 3)

### ill-fitting normal model
nd$y <- rnorm(nrow(nd), mean = c(0, .5)[nd$groups], sd = c(1, 1.5)[nd$groups])
with(with(nd, split(y, groups)),
     ppplot(G1, G2, conf.level = .95,
            conf.args = list(link = "probit", type = "Wald", col = 2)))
### inappropriate probit model
coef(free1way(y ~ groups, data = nd, link = "probit"))

}
\keyword{hplot}
\keyword{distribution}
