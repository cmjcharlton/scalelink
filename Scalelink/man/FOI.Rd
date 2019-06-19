\docType{data}
\name{FOI}
\alias{FOI}
\title{File of interest}
\format{A data frame with 7742 observations on the following 5 variables:
  \describe{
    \item{\code{id}}{Record Identifier (not used for linking).}
    \item{\code{Day}}{Day of Week.}
    \item{\code{Month}}{Month of Year.}
    \item{\code{Surname}}{Surname.}
    \item{\code{Year}}{Year.}
    \item{\code{Sex}}{Gender: with codes \code{1} Male and \code{2} Female.}
}}
\source{
Synthetic data created by Harvey Goldstein
}

\description{
  File of interest data data with 7742 records and 5 variables.}
\details{
  The \code{FOI} dataset is one of the sample datasets provided with this package for demonstration purposes.}

\examples{
data(FOI, package = "Scalelink")
summary(FOI)
}
\keyword{datasets}
