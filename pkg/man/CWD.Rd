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

Coarse woody debris (CWD, dead wood greater 10 cm diameter) is a large stock
of carbon in tropical forests, yet the flux of carbon out of this pool, via
respiration, is poorly resolved (Chambers et al., 2001).
The heterotrophic process involved in CWD respiration should respond to
reductions in moisture availability, which occurs during dry season
(Chambers et al., 2001).  CWD respiration measurements
were taken in a tropical forest in west French Guiana, which experiences
extreme contrasts in wet and dry season
(Bonal et al., 2008).  An infrared gas analyzer and a clear
chamber sealed to the wood surface were used to measure the flux of carbon
out of the wood (Stahl et al., 2001).
Measurements were repeated on six pieces of wood 13 times from July to November
2011, during the transition into and out of the dry season.  The aim is to
assess if there were shifts in the CWD respiration of any of the pieces
in response to the transition into (early August) and
out of (late October) the dry season. 


}
\references{

Jeffrey Q. Chambers, Joshua P. Schimel, Antonio D. Nobre (2001).
Respiration from Coarse Wood Litter in Central Amazon Forests,
\emph{Biogeochemistry},\bold{52}(2), 115--131.

Damien Bonal, Alexandre Bosc, St\'{e}phane Ponton,
Jean-Yves Goret, Beno\^{i}t Burban, Patrick Gross,
Jean-Marc Bonnefond, Jan Elbers, Bernard Longdoz, 
Daniel Epron, Jean-Marc Guehl, Andr\'{e} Granier (2008).
Impact of Severe Dry Season on Net Ecosystem Exchange
in the Neotropical Rainforest of French Guiana.
\emph{Global Change Biology}, \bold{14}(8), 1917--1933.

Cl\'{e}ment Stahl and Beno\^{i}t Burban and Jean-Yves Goret
and Damien Bonal (2011). Seasonal Variations in Stem {CO$_2$} Efflux in the
Neotropical Rainforest of French Guiana.
\emph{Annals of Forest Science}, \bold{68}(4), 771--782.

}

\source{
The coarse woody debris respiration data were kindly provided by Lucy Rowland (School of
GeoSciences, University of Edinburgh).
}
\examples{

  data("CWD")
  CWD[1:6] <- 100 * CWD[1:6] # scaling (to avoid warning)
  maxstat_test(sample2 + sample3 + sample4 +
               sample6 + sample7 + sample8 ~ trend, data = CWD,
               distribution = approximate(1e5))

}
\keyword{datasets}
