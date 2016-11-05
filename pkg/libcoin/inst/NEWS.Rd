
\name{NEWS}
\title{NEWS file for the libcoin package}

\section{Changes in Version 0.6-0 (2016-09-01)}{
  \itemize{
    \item{Give meaningful defaults and avoid calling \code{missing()} in
      \code{LinStatExpCov()}.}
    \item{Add \code{call_libcoin()}. This is a dirty by-pass and must not be
      used under any circumstances.}
  }
}

\section{Changes in Version 0.5-1 (2016-08-23)}{
  \itemize{
    \item{Export \code{ctabs()} function for efficient cross-tabulation.}
  }
}

\section{Changes in Version 0.5-0 (2016-08-17)}{
  \itemize{
    \item{\pkg{libcoin} passes \code{R CMD check} and reproduces all \pkg{coin}
      results.}
  }
}
