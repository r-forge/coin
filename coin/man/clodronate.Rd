\name{clodronate}
\alias{clodronate}
\docType{data}
\title{ Clodronate and Multiple Myeloma }
\description{
  Prospective Randomized Trial of Clodronate in Patients with Multiple
  Myeloma
}
\usage{data(clodronate)}
\format{
  A data frame with 59 observations on the following 3 variables.
  \describe{
    \item{BoneRes6}{Bone response after 6 month. A factor with levels 
                    \code{NC} (no change), \code{PD} (progressing)
                    and \code{R} (receding). }
    \item{BoneRes12}{Bone response after 12 month.}
    \item{group}{a factor with levels \code{control} (control group) 
                 and \code{treat} (treatment group receiving clodronate).}}
}
\details{
  
}
\source{
  Data kindly provided by Dr. L. Edler. 
}
\references{
  M.E. Heim, M.R. Clemens, W. Queisser, M. Pecherstorfer, Ch. Boewer, 
  M. Herold, A. Franke, Z. Herrmann, R. Loose, L. Edler (1995). Prospective
  randomized trial of dichloromethylene bisphosphonate (Clodronate) in
  patients with multiple myeloma requiring treatment. A multicenter study.
  \emph{Onkologie}, \bold{18}, 439--448.
}
\examples{

### <FIXME> Linear-by-Linear association test!!!
data(clodronate)
cmh_test(BoneRes6 ~ group, data = clodronate, distribution = "approx")
### <FIXME>

}
\keyword{datasets}
