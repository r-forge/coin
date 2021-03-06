\name{coin}
\alias{coin}
\title{General Information on the coin Package}
\description{

The \pkg{coin} package implements a general framework for conditional
inference procedures, commonly known as \emph{permutation tests},
theoretically derived by Strasser & Weber (1999). The conditional
expectation and covariance for a broad class of multivariate linear
statistics as well as the corresponding multivariate limiting distribution
was derived by Strasser & Weber (1999). These results are 
utilized to construct tests for independence between two sets of variables. 

Beside a general implementation of the abstract framework the package offers
a rather huge set of convenience functions implementing well known classical
as well as less prominent classical and non-classical test procedures 
in a conditional inference framework. Examples are linear rank statistics for
the two- and K-sample location and scale problem against ordered and
unordered alternatives including post-hoc tests for arbitrary contrasts,
tests of independence for contingency tables, two- and K-sample tests for
censored data, tests for independence of two continuous variables as well as 
tests for marginal homogeneity and symmetry. Conditional counterparts of 
most of the classical procedures
given in famous text books like Hollander & Wolfe (1999) or Agresti (2002)
can be implemented as part of the general framework without much effort.
Approximations of the exact null distribution via the limiting distribution and 
conditional Monte-Carlo procedures are available for every test while the 
exact null distribution is currently available for two-sample problems only. 

}
\author{

    This package is written by \cr
    Torsten Hothorn <Torsten.Hothorn@R-project.org>, \cr 
    Kurt Hornik <Kurt.Hornik@R-project.org>, \cr
    Mark van de Wiel <m.a.v.d.wiel@TUE.nl>, \cr
    Henric Winell <henric.winell@icm.uu.se> and \cr
    Achim Zeileis <Achim.Zeileis@R-project.org>.

}
\references{

    Helmut Strasser & Christian Weber (1999).
    On the asymptotic theory of permutation statistics,
    \emph{Mathematical Methods of Statistics} \bold{8}, 220--250.  

    Myles Hollander & Douglas A. Wolfe (1999).
    \emph{Nonparametric Statistical Methods, 2nd Edition}.
    New York: John Wiley & Sons.

    Alan Agresti (2002). \emph{Categorical Data Analysis} Hoboken, New
    Jersey: John Wiley & Sons.

    Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel & Achim Zeileis (2006).
    A Lego system for conditional inference, \emph{The American Statistician},
    \bold{60}(3), 257--263.

    Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel & Achim Zeileis (2008).
    Implementing a class of permutation tests: The coin package,
    \emph{Journal of Statistical Software}, \bold{28}(8), 1--23.
    \url{http://www.jstatsoft.org/v28/i08/}

}
\examples{

    ### if you are interested in the internals:
    ### generate doxygen documentation
    \dontrun{

        ### download src package into temp dir
        tmpdir <- tempdir()
        tgz <- download.packages("coin", destdir = tmpdir)[2]
        ### extract
        untar(tgz, exdir = tmpdir)
        wd <- setwd(file.path(tmpdir, "coin"))
        ### run doxygen (assuming it is there)
        system("doxygen doxygen.cfg")
        setwd(wd)
        ### have fun
        browseURL(file.path(tmpdir, "coin", "inst",
                            "documentation", "html", "index.html"))
    }
}
\keyword{misc}
