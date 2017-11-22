\name{cenwave}
\alias{cenwave}
\docType{data}
\title{
Filter Central Wavelengths
}
\description{
Table of central wavelengths of all provided filters.
}
\usage{data("cenwave")}
\format{
  A data frame with 26 observations on the following 2 variables.
  \describe{
    \item{\code{filter}}{A character vector}
    \item{\code{cenwave}}{A numeric vector}
  }
}
\examples{
data('cenwave')
print(cenwave)

barplot(height=cenwave$cenwave, names.arg=cenwave$filter, log='y')
}
\keyword{datasets}
