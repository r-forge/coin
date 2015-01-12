\name{GTSG}
\alias{GTSG}
\docType{data}
\title{Gastrointestinal Tumor Study Group}

\encoding{UTF-8}

\description{
  A randomized clinical trial in gastric cancer.
}
\usage{data("GTSG")}
\format{
  A data frame with observations from the GTSG study.
  \describe{
    \item{time}{a numeric, survival time in days.}
    \item{event}{a numeric, status indicator: 0 for censored
      observations and 1 otherwise.}
    \item{group}{a factor with levels \code{"Chemotherapoy+Radiation"}
      and \code{"Chemotherapy"}.}
  }
}
\details{
  A clinical trial comparing chemotherapy alone versus a combination of
  chemotherapy and radiation therapy in the treatment of locally
  advanced, nonresectable gastric carcinoma.
}
\note{
  There is substantial separation between the estimated survival
  distributions at 8 to 10 months, but by month 26 the distributions
  intersect.
}
\source{
  Donald M. Stablein, Walter H. Carter, Jr., & Joel W. Novak (1981).
  Analysis of Survival Data with Nonproportional Hazard Functions.
  \emph{Controlled Clinical Trials} \bold{2}(2), 149--159.
}
\references{
  Thierry Moreau, Jean Maccario, Joseph Lellouch & Catherine Huber
  (1992).
  Weighted Log Rank Statistics for Comparing Two Distributions.
  \emph{Biometrika} \bold{79}(1), 195--198.

  Pascale Tubert-Bitter, Andrew Kramar, Jean Jacques \enc{Chalé}{Chale}
  & Thierry Moureau (1994).
  Linear Rank Tests for Comparing Survival in Two Groups with Crossing
  Hazards.
  \emph{Computational Statistics & Data Analysis} \bold{18}(5),
  547--559.

  Wei Shen & Chap T. Le (2000).
  Linear Rank Tests for Censored Survival Data.
  \emph{Communications in Statistics -- Simulation and Computation}
  \bold{29}(1), 21--36.
}
\examples{
  ### plot Kaplan-Meier estimates
  plot(survfit(Surv(time / (365.25 / 12), event) ~ group, data = GTSG),
       lty = 1:2,
       ylab = "\% Survival", xlab = "Survival Time in Months")
  legend("topright", lty = 1:2,
         c("Chemotherapy+Radiation", "Chemotherapy"), bty = "n")

  ### logrank test
  surv_test(Surv(time, event) ~ group, data = GTSG)

  ### Prentice test
  surv_test(Surv(time, event) ~ group, data = GTSG,
            type = "Prentice")

  ### testing against Weibull-type alternatives, see Moreau et al (1992)
  moreau_weight <- function(time, n.risk, n.event)
      1 + log(-log(cumprod(n.risk / (n.risk + n.event))))

  independence_test(Surv(time, event) ~ group, data = GTSG,
                    ytrafo = function(data)
                        trafo(data, surv_trafo = function(y)
                            logrank_trafo(y, weight = moreau_weight)))

  ### testing against crossing-curve alternatives, see Shen & Le (2000)
  shen_trafo <- function(x)
      ansari_trafo(logrank_trafo(x, type = "Prentice"))

  independence_test(Surv(time, event) ~ group, data = GTSG,
                    ytrafo = function(data)
                        trafo(data, surv_trafo = shen_trafo))
}
\keyword{datasets}
