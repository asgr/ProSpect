\name{photom}
\alias{photom}
\alias{Lum2Flux}
\alias{Lum2FluxFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral Processing for Photometry
}
\description{
These functions can manipulate restframe spectra in useful ways, e.g. placing it at a target redshift with the correct observed luminosity and band stretching.
}
\usage{
photom(wave, lum, filters = "all", z = 0.1, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = "planck")
Lum2Flux(wave, lum, z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR,
OmegaR = 0, w0 = -1, wprime = 0, ref='planck')
Lum2FluxFactor(z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR,
OmegaR = 0, w0 = -1, wprime = 0, ref='planck')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; the spectral wavelength in Angstroms. For convenience, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{lum}.
}
  \item{lum}{
Numeric vector; the spectral luminosity in Lsun / Angstrom (BC03lr / BC03hr and EMILES format).
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters.
}
  \item{z}{
Numeric scalar; cosmological redshift, where z must be > -1.
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.3).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
}
  \item{OmegaR}{
Numeric scalar; Omega Radiation today (default is 0, but OmegaM/3400 is typical).
}
  \item{w0}{
The value of dark energy equation of state at z=0. See \code{\link{cosgrow}} for more details.
}
  \item{wprime}{
The evolution term that governs how the dark energy equation of state evolves with redshift. See \code{\link{cosgrow}} for more details.
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck=Planck15 and WMAP=WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
}
}
\details{
\code{photom} is a useful utility function since it fully converts between an intrinsic spectrum (e.g. of the type created when integrating a star formation history over differ stellar population models) all the way to a band-pass magnitude in AB mags. The latter is how many surveys represent their photometry, so this is usually a sensible unit space to use for feature fitting etc. In some circumstances fluxes (CGS or Jansky type) might be more appropriate, espeically when signal to noise is low a meaurements might be negative.
}
\value{
\item{photom}{Data frame; first column is the filter name, second columns is the central wavelength, third column is the AB mag.}
\item{Lum2Flux}{Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux (erg/s/cm^2/Angstrom).}
\item{Lum2FluxFactor}{Numeric scalar; the scaling factor to convert between intrinsic Lsun/Angstrom and erg/s/cm^2/Angstrom.}
}
\references{
Bruzual and Charlot, 2003, MNRAS, 344, 1000
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFHp4}}
}
\examples{
#Make a 10^10 solar mass SSP that is 5 Grys old and solar metallicity:

#First with BC03lr:
BC03lr$Age[161]

photoBC03=photom(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
plot(photoBC03[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)

#Next with EMILES
EMILES$Age[35]

photoEMILES=photom(EMILES$Wave, EMILES$Zspec[[9]][35,]*1e10)
plot(photoEMILES[,2:3], log='x', ylim=c(26,17), xlab=EMILES$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)

#We can show the flux directly too:

fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ magnitude }% use one of  RShowDoc("KEYWORDS")
\keyword{ photometry }% __ONLY ONE__ keyword per line
