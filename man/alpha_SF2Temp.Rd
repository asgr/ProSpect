\name{alpha_SF2Temp}
\alias{alpha_SF2Temp}
\alias{alpha_SF2Temp_func}
\alias{Temp2alpha_SF_func}
\docType{data}
\title{
Dale alpha_SF to Temperature Conversion Table and Functions
}
\description{
This is a by eye Dale alpha_SF to temperature (K) conversion table (\code{alpha_SF2Temp}) and associated functions ((\code{alpha_SF2Temp_func} / \code{Temp2alpha_SF_func})) that let you convert in both directions. Here we fit to the Dale templates with no AGN contribution, and use a grey body for the fit.
}
\usage{
data("alpha_SF2Temp")
data("alpha_SF2Temp_func")
data("Temp2alpha_SF_func")
}
\format{
  A data frame with 7 observations on the following 3 variables.
  \describe{
    \item{\code{alpha_SF}}{A numeric vector; alpha slope of the radiation field in a Dale dust template.}
    \item{\code{Temp}}{A numeric vector; approximate grey body temperature of the Dale dust template.}
    \item{\code{Beta}}{A numeric vector; approximate grey body beta of the Dale dust template.}
  }
}
\details{
Useful for a rough temperature (easier to think of the dust in these terms. In detail this will not be correct, since the dust is heated by a continuum, and there is no single appropriate temperature to state. However, adjusting also the grey body beta gives us a good rough fit that is indicative at least.
}
\references{
Dale et al, 2014, ApJ, 784, 11
}
\examples{
data(alpha_SF2Temp)
data(alpha_SF2Temp_func)
plot(alpha_SF2Temp[,'alpha_SF'], alpha_SF2Temp[,'Temp'], xlab='alpha_SF', ylab='Temp/K')
curve(alpha_SF2Temp_func, 0, 4, add=TRUE)
}
\keyword{datasets}
