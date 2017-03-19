\docType{data}
\name{LDF}
\alias{LDF}
\title{Linking data file}
\format{A data frame with 10000 observations on the following 6 variables:
  \describe{
    \item{\code{id}}{Identifier.}
    \item{\code{Day}}{Day of Week.}
    \item{\code{Month}}{Month of Year.}
    \item{\code{Year}}{Year.}
    \item{\code{Val}}{Value.}
    \item{\code{Sex}}{Gender: with codes \code{1} Male and \code{2} Female.}
}}
\source{
Synthetic data created by Harvey Goldstein
}

\description{
  Linking data file data with 10000 records and 5 variables.}
\details{
  The \code{LDF} dataset is one of the sample datasets provided with this package for demonstration purposes.
  This version include records with missing data
}

\examples{
\dontrun{

data(LDF, package = "Scalelink")
summary(LDF)
}
}
\keyword{datasets}
