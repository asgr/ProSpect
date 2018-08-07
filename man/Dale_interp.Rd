\name{Dale_interp}
\alias{Dale_interp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate the Dale 2014 Templates
}
\description{
Interpolates between non-computed values of \option{AGNfrac} and \option{alpha_SF}. This is useful when creating a generic generative model of a galaxy SED, where it can harm the sampling to have strictly discretised values for variables we wish to sample over.
}
\usage{
Dale_interp(alpha_SF = 1.5, AGNfrac = 0, type = "Msol", Dale = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha_SF}{
Numeric scalar; desired interpolated alpha slope of the star forming population.
}
  \item{AGNfrac}{
Numeric scalar; desired interpolated AGN fraction linearly mixed in energy over the range 5-20 microns (see Dale 2014).
}
  \item{type}{
Character scalar; type of Dale template to interpolate. One of: Orig, Msol, NormTot, NormAGN, NormSFR (see \code{\link{Dale}}).
}
  \item{Dale}{
Pass in the Dale dust library directly. Must be one of \code{\link{Dale_Orig}}, \code{\link{Dale_Msol}}, \code{\link{Dale_NormTot}}, \code{\link{Dale_NormAGN}}, \code{\link{Dale_NormSFR}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
}
\details{
This function does a simple 2D bilinear interpolation to create templates for values of \option{AGNfrac} and \option{alpha_SF} that are not precisely provided by \code{\link{Dale}}.
}
\value{
Two column data.frame containing the wavelength (Wave) and the interpolated luminosity values (Aspec) in the same units as the input Dale dust library (from argument \option{Dale}).
}
\references{
Dale et al, 2014, ApJ, 784, 11
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Dale}}, \code{\link{interp_param}}
}
\examples{
#Example output showing also the 4 surrounding interpolation points used.
temp=Dale_interp(1.53125, 0.025)

plot(temp, log='x', type='l', col='red')
lines(Dale_Orig$Wave, Dale_Msol$Aspec[[1]][24,])
lines(Dale_Orig$Wave, Dale_Msol$Aspec[[1]][25,])
lines(Dale_Orig$Wave, Dale_Msol$Aspec[[2]][24,])
lines(Dale_Orig$Wave, Dale_Msol$Aspec[[2]][25,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Dale }% use one of  RShowDoc("KEYWORDS")
\keyword{ Dust }% __ONLY ONE__ keyword per line
