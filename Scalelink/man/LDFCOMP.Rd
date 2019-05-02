\docType{data}
\name{LDFCOMP}
\alias{LDFCOMP}
\title{Linking data file}
\format{A data frame with 8142 observations on the following 5 variables:
  \describe{
    \item{\code{id}}{Record Identifier (not used for linking).}
    \item{\code{Day}}{Day of Week.}
    \item{\code{Month}}{Month of Year.}
    \item{\code{Year}}{Year.}
    \item{\code{Sex}}{Gender: with codes \code{1} Male and \code{2} Female.}
}}
\source{
Synthetic data created by Harvey Goldstein
}

\description{
  File of interest data data with 8142 records and 5 variables.}
\details{
  The \code{LDFCOMP} dataset is one of the sample datasets provided with this package for demonstration purposes.
  This version has records containing missing data removed
}

\examples{
data(LDFCOMP, package = "Scalelink")
summary(LDFCOMP)
}
\keyword{datasets}
