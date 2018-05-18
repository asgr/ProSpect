\name{magAB2Jansky}
\alias{magAB2Jansky}
\alias{Jansky2magAB}
\alias{magAB2CGS}
\alias{CGS2magAB}
\alias{Jansky2CGS}
\alias{CGS2Jansky}
\alias{AbsoluteToWHz}
\alias{WHzToAbsolute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between AB mag, Jansky and CGS
}
\description{
Simple functions to convert between the common flux and mag systems.
}
\usage{
magAB2Jansky(x)
Jansky2magAB(x)
magAB2CGS(x)
CGS2magAB(x)
Jansky2CGS(x)
CGS2Jansky(x)
AbsoluteToWHz(x)
WHzToAbsolute(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the value to be converted (in units of the LHS of the function name).
}
}
\details{
For clarity the conversion used are provided below:

magAB2Jansky = 10^(-0.4*(magAB-8.9))

Jansky2magAB = -2.5*log10(Jansky)+8.9

magAB2CGS = 10^(-0.4*(magAB+48.6))

CGS2magAB = -2.5*log10(CGS)-48.6

Jansky2CGS = Jansky*1e-23

CGS2Jansky = CGS*1e23

Explictly for absolute magnitudes we can convert to the W/Hz convention. This cannot be used for general magnitudes (or at least, the outputs will not make much sense).

AbsoluteToWHz = 4*pi*((3.086e+17)^2)*1e-26*magAB2Jansky(Absolute)

WHzToAbsolute = Jansky2magAB(WHz/(4*pi*((3.086e+17)^2)*1e-26))

}
\value{
Numeric vector; the converted value (in units of the RHS of the function name).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{magABcalc}}
}
\examples{
magAB2Jansky(20)
Jansky2magAB(magAB2Jansky(20))
Jansky2CGS(magAB2Jansky(20))

magAB2CGS(20)
CGS2magAB(magAB2CGS(20))
CGS2Jansky(magAB2CGS(20))

AbsoluteToWHz(-22)
WHzToAbsolute(AbsoluteToWHz(-22))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ magAB }% use one of  RShowDoc("KEYWORDS")
\keyword{ Jansky }% __ONLY ONE__ keyword per line
\keyword{ CGS }
