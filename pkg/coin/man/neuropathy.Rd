\name{neuropathy}
\alias{neuropathy}
\docType{data}
\title{ Acute Painful Diabetic Neuropathy }
\description{
  The logarithm of the ratio of pain scores at baseline and after four
  weeks for a control and treatment group.
}
\usage{data("neuropathy")}
\format{
  A data frame with 58 observations on the following 2 variables.
  \describe{
    \item{pain}{pain scores: ln(baseline/final).}
    \item{group}{a factor with levels \code{control} and \code{treat}.}
  }
}
\details{
  Data from Table 1 of Conover & Salsburg (1988).
}
\source{

  William J. Conover & David S. Salsburg (1988).
  Locally Most Powerful Tests for Detecting Treatment Effects When Only
  a Subset of Patients Can Be Expected to \dQuote{Respond} to Treatment.
  \emph{Biometrics} \bold{44}(1), 189--196.

}
\examples{

  ### compare with Table 2 of Conover & Salsburg (1988)
  oneway_test(pain ~ group, data = neuropathy, alternative = "less",
              distribution = "exact")

  wilcox_test(pain ~ group, data = neuropathy, alternative = "less",
              distribution = "exact")

  oneway_test(pain ~ group, data = neuropathy, alternative = "less",
              distribution = approximate(B = 10000),
              ytrafo = function(data)
                  trafo(data, numeric_trafo = consal_trafo))

  ### maximum-type test for a range of 'a' values
  it <- independence_test(pain ~ group, data = neuropathy,
                          alternative = "less",
                          distribution = approximate(B = 10000),
                          ytrafo = function(data)
                              trafo(data, numeric_trafo = function(y)
                                  consal_trafo(y, a = 2:7)))
  pvalue(it, method = "single-step")

}
\keyword{datasets}
