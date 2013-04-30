\name{CWD}
\alias{CWD}
\docType{data}
\title{
Coarse Woody Debris Data
}
\description{
Carbon flux on six pieces of wood.
}
\usage{data(CWD)}
\format{
  A data frame with 13 observations on the following 8 variables.
  \describe{
    \item{\code{sample2}}{carbon flux measurement for 2nd piece of wood.}
    \item{\code{sample3}}{carbon flux measurement for 3rd piece of wood.}
    \item{\code{sample4}}{carbon flux measurement for 4th piece of wood.}
    \item{\code{sample6}}{carbon flux measurement for 6th piece of wood.}
    \item{\code{sample7}}{carbon flux measurement for 7th piece of wood.}
    \item{\code{sample8}}{carbon flux measurement for 8th piece of wood.}
    \item{\code{trend}}{measurement day (in days from beginning).}
    \item{\code{time}}{date of measurement.}
  }
}
\details{

Aninfra red gas analyser and a clear chamber sealed to the wood surface were
used to measure the flux of carbon out of the wood. 
Measurements were repeated on CWD samples 13 times from July to November
2011, during the transition into and out of the dry season.  The aim was to
assess if there was shifts in the CWD respiration in response to the shift
in the seasons.
}
\source{

Cl\'{e}ment Stahl and Beno\^{i}t Burban and Jean-Yves Goret
and Damien Bonal (2011). Seasonal Variations in Stem {CO$_2$} Efflux in the
Neotropical Rainforest of French Guiana.
\emph{Annals of Forest Science}, \bold{68}(4), 771--782.
}
\examples{

  data("CWD")
  CWD[1:6] <- 100 * CWD[1:6] # scaling (to avoid warning)
  maxstat_test(sample2 + sample3 + sample4 +
               sample6 + sample7 + sample8 ~ trend, data = CWD,
               distribution = approximate(1e5))

}
\keyword{datasets}
