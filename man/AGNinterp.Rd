\name{AGNinterp}
\alias{AGNinterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate the Fritz 2006 AGN Templates
}
\description{
Interpolates between non-computed values of the Fritz 2006 AGN Templates using a mixture of schemes. This is similar to the approach used for \code{\link{Dale_interp}}. The model is much too complex to summarise here, so interested readers should refer to the source papers if using this function for science.
}
\usage{
AGNinterp(lum = 1e+44, ct = 60, al = 4, be = -0.5, ta = 1, rm = 60, an = 30,
Fritz = NULL)
}
\arguments{
  \item{lum}{
Numeric scalar; AGN bolometric luminosity in erg/s. This is of the illuminating accretion source, so actual observed output will vary depending on the geometry of the dust torus.
}
  \item{ct}{
Numeric scalar; opening angle of torus in degrees.
}
  \item{al}{
Numeric scalar; gamma dust parameter.
}
  \item{be}{
Numeric scalar; bBeta dust parameter.
}
  \item{ta}{
Numeric scalar; optical depth tau.
}
  \item{rm}{
Numeric scalar; outer to inner torus radius ratio.
}
  \item{an}{
Numeric scalar; angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned.
}
  \item{Fritz}{
The \option{Fritz} data. Providing this directly reduces disk IO since it will use the copy already loaded into RAM. If not provided it will load the data from disk.
}
}
\details{
The below is taken directly from the included documentation:

A composition of power-laws, with different spectral indexes as a function of the wavelength, was used to model the primary source (i.e. accretion disk SED) SED taken from the file pwlbb05.dat.                                          

Standard, Galactic composition of dust, with graphite and astronomical silicates. The grain size distribution is a power law with index -3.5. Scattering and absorption coefficients taken from Laor and Draine (1993). Number of different species: 2 (Graphite and Silicate only).                      

The inner radius is set by the sublimation temperature of graphite grains: T = 1500 K. The dust sublimation radius is Rmin = 1.3*L46^0.5*T^(-2.8) pc [Barvainis (1987)] where L46=L/(1e46 erg/s), and T is the sublimation temperature. The dust density throughout the torus is expressed by the following: rho(r,theta)=A*r^(beta)exp(-gamma*cos(theta)), where A is a normalization costant, beta and gamma are, in the name of the model, be and al, respectively.
}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed luminosity (erg/s/Ang).
}
\references{
Feltre et al, 2012, MNRAS, 426, 120 \cr
Fritz et al, 2006, MNRAS, 366, 767
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Fritz}}
}
\examples{
plot(AGNinterp(), type='l', log='xy', xlab='Wave / Ang', ylab='Lum / erg/s/A')
lines(AGNinterp(an=0), col='blue')
lines(AGNinterp(an=90), col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Fritz }% use one of  RShowDoc("KEYWORDS")
