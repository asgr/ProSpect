\name{Dale_M2L}
\alias{Dale_M2L}
\alias{Dale_M2L_func}
\docType{data}
\title{
Dale Template Mass-to-Light
}
\description{
Returns the correct mass-to-light for a given \option{alpha_SF} (which is akin to dust temperature).
}
\usage{
data("Dale_M2L")
data("Dale_M2L_func")
}
\format{
  \code{Dale_M2L} is a data frame with 64 observations on the following 2 variables.
  \describe{
    \item{\code{alpha_SF}}{Numeric vector; desired interpolated alpha slope of the star forming population. Lower values mean hotter dust.}
    \item{\code{M2L}}{Numeric vector; mass-to-light of the dust.}
  }
  \code{Dale_M2L_func} is a spline fit to \code{Dale_M2L}, and takes the input \option{alpha_SF} and returns the respective \option{M2L}.
}
\details{
Simple data set and corresponding function that provide the mass-to-light which is used in \code{\link{ProSpectSED}} to compute dust masses.
}
\examples{
data(Dale_M2L)
data(Dale_M2L_func)
plot(Dale_M2L, log='y')
curve(Dale_M2L_func, col='red', add=TRUE)

data(Dale_Msol)
sum(c(0,diff(Dale_Msol$Wave))*Dale_Msol$Aspec[[1]][30,])
Dale_M2L[30,'M2L']
}
\keyword{datasets}
