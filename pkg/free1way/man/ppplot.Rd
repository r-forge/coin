\name{ppplot}
\alias{ppplot}
\title{
Probability-probability Plots
}
\description{
\code{ppplot} produces a PP plot of two numerical variables, essentially
the Receiver Operating Characteristic (ROC) curve. 
If \code{conf.level} is given, a confidence band for the ROC curve under a
distribution-free model is plotted.
}
\usage{
ppplot(x, y, plot.it = TRUE, xlab = deparse1(substitute(x)), 
       ylab = deparse1(substitute(y)), ..., conf.level = NULL, 
       conf.args = list(link = "logit", type = "Wald", col = NA, border = NULL))
}
\arguments{
  \item{x}{The first sample for \code{ppplot}, a numerical variable.
}
  \item{y}{The second sample for \code{ppplot}, a numerical variable.
}
  \item{plot.it}{logical. Should the result be plotted?
}
  \item{xlab, ylab}{The \code{xlab} and \code{ylab} refer to the y
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

### well fitting logistic model
nd <- data.frame(groups = gl(2, 50, labels = paste0("G", 1:2)))
nd$y <- rlogis(nrow(nd), location = c(0, 2)[nd$groups])
with(with(nd, split(y, groups)),
     ppplot(G1, G2, conf.level = .95,
            conf.args = list(link = "logit", type = "Wald", col = 2)))
### with appropriate Wilcoxon test and log-odds ratio
coef(free1way(y ~ groups, data = nd))

### ill fitting normal model
nd$y <- rnorm(nrow(nd), mean = c(0, .5)[nd$groups], sd = c(1, 1.5)[nd$groups])
with(with(nd, split(y, groups)),
     ppplot(G1, G2, conf.level = .95,
            conf.args = list(link = "probit", type = "Wald", col = 2)))
### inappropriate probit model
coef(free1way(y ~ groups, data = nd, link = "probit"))

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
