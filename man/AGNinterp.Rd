\name{AGNinterp}
\alias{AGNinterp}
\alias{Fritz_interp}
\alias{SKIRTOR_interp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate the Fritz (2006) and SKIRTOR (2016) AGN Templates
}
\description{
Interpolates between non-computed values of the Fritz 2006 AGN Templates using a mixture of schemes. This is similar to the approach used for \code{\link{Dale_interp}}. The model is much too complex to summarise here, so interested readers should refer to the source papers if using this function for science.
}
\usage{
Fritz_interp(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, al = 4, be = -0.5,
  Fritz = NULL)

SKIRTOR_interp(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
  SKIRTOR = NULL)
}
\arguments{
  \item{lum}{
Numeric scalar; AGN bolometric luminosity in erg/s. This is of the illuminating accretion source, so actual observed output will vary depending on the geometry of the dust torus.
}
  \item{ct}{
Numeric scalar; opening angle of torus in degrees.
}
  \item{rm}{
Numeric scalar; outer to inner torus radius ratio.
}
  \item{an}{
Numeric scalar; angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned. See Notes.
}
  \item{ta}{
Numeric scalar; optical depth tau.
}
  \item{al}{
Numeric scalar; gamma dust parameter.
}
  \item{be}{
Numeric scalar; beta dust parameter.
}
  \item{p}{
Numeric scalar; Index for radial density gradient.
}
  \item{q}{
Numeric scalar; Index for angular density gradient.
}
  \item{Fritz}{
The \option{Fritz} data. Providing this directly reduces disk IO since it will use the copy already loaded into RAM. If not provided it will load the data from disk.
}
  \item{SKIRTOR}{
The \option{SKIRTOR} data. Providing this directly reduces disk IO since it will use the copy already loaded into RAM. If not provided it will load the data from disk.
}
}
\details{
The below is taken directly from the included Fritz AGN documentation:

A composition of power-laws, with different spectral indexes as a function of the wavelength, was used to model the primary source (i.e. accretion disk SED) SED taken from the file pwlbb05.dat.

Standard, Galactic composition of dust, with graphite and astronomical silicates. The grain size distribution is a power law with index -3.5. Scattering and absorption coefficients taken from Laor and Draine (1993). Number of different species: 2 (Graphite and Silicate only).

The inner radius is set by the sublimation temperature of graphite grains: T = 1500 K. The dust sublimation radius is Rmin = 1.3*L46^0.5*T^(-2.8) pc [Barvainis (1987)] where L46=L/(1e46 erg/s), and T is the sublimation temperature. The dust density throughout the torus is expressed by the following: rho(r,theta)=A*r^(beta)exp(-gamma*cos(theta)), where A is a normalization constant, beta and gamma are, in the name of the model, be and al, respectively.

The below is taken directly from the included SKIRTOR AGNdocumentation:

We modeled the distribution of the dust in the torus as a two-phase medium, consisting of a large number of high-density clumps embedded in a smooth dusty component of low density.

Previously, we have found that an advantageous property of such two-phase medium is that it can produce attenuated silicate features, while at the same time have a pronounced NIR emission, which is challenging for both smooth and clumpy models alone.
}
\note{
If \option{an} this is larger then \option{ct}/2 then by definition a significant component of the AGN accretion disk light is being attenuated by the torus (AKA type II AGN). If \option{an} is less then \option{ct}/2 then we have a direct line-of-sight to the AGN accretion disk and will observe much more UV flux in particular (AKA type I AGN, including blazars and quasars).
}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed luminosity (erg/s/Ang).
}
\references{
Feltre et al, 2012, MNRAS, 426, 120 \cr
Fritz et al, 2006, MNRAS, 366, 767 \cr
Stalevski et al, 2012, MNRAS, 420, 2756 \cr
Stalevski et al, 2016, MNRAS, 458, 2288
}
\author{
Aaron Robotham & Adam Marshall
}
\seealso{
\code{\link{Fritz}}, \code{\link{SKIRTOR}}
}
\examples{
\dontrun{
library(magicaxis)
colvec = hcl.colors(10, 'hawaii')

magplot(NA, NA, log='xy', xlab='Wave / Ang', ylab='Lum / erg/s/A',
        xlim=c(1e2,1e7), ylim=c(1e30,1e42))
for(i in 1:10){
  lines(Fritz_interp(an=(i-1)*10), col=colvec[i], lty=2, lwd=2)
  lines(SKIRTOR_interp(an=(i-1)*10), col=colvec[i], lwd=2)
}

magbar('topright', range=c(0,90), col=colvec, labN=3, prettybase=30, title='an (Deg)')
legend('bottomleft', legend=c('SKIRTOR', 'Fritz'), lty=c(1,2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Fritz }% use one of  RShowDoc("KEYWORDS")
\keyword{ SKIRTOR }% use one of  RShowDoc("KEYWORDS")
