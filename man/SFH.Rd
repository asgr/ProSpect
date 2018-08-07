\name{Star Formation}
\alias{SFH}
\alias{SFHp4}
\alias{SFHp5}
\alias{SFHfunc}
\alias{SM}
\alias{SMstarp4}
\alias{SMstarp5}
\alias{SMstarfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Various Star Formation History Functions
}
\description{
Computes output galaxy spectra  and stellar masses using variable star formation histories. This can be either distinct phases with \code{SFHp4} and \code{SFHp5}, or a functional form with \code{SFHfunc}. Uses either BC03 or EMILES libraries. For \code{SFHp4} and \code{SFHp5}, each star formation phase can have a different metallicity associated with it (default solar).
}
\usage{
SFHp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop = "BC03lr", speclib = NULL, tau_birth = 1,
tau_screen = 0.3, filters = "all", Z = c(5, 5, 5, 5), z = 0.1, H0 = 67.8, OmegaM = 0.308,
OmegaL = 1 - OmegaM, outtype = "mag", cossplit = c(9e+09, 1.3e+10), dosplit = FALSE,
sparse = 1)

SFHp5(burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10, oldmass = 1e+10,
ancientmass = 1e+10, burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09),
midage = c(1e+09, 5e+09), oldage = c(5e+09, 9e+09), ancientage = c(9e+09, 1.3e+10),
stellpop = "BC03lr", speclib = NULL, tau_birth = 1, tau_screen = 0.3, filters = "all",
Z = c(5, 5, 5, 5, 5), z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM,
outtype = "mag", cossplit = c(9e+09, 1.3e+10), dosplit = FALSE, sparse = 1)

SFHfunc(massfunc = function(age, SFR = 1) { ifelse(age < 1e+10, SFR, 0)},
forcemass = FALSE, unimax = 1.38e+10, agescale = 1, stellpop = "BC03lr", speclib = NULL,
tau_birth = 1, tau_screen = 0.3, filters = "all", Z = 5, z = 0.1, H0 = 67.8,
OmegaM = 0.308, OmegaL = 1 - OmegaM, outtype = "mag", sparse = 1, intSFR=FALSE, ...) 

SMstarp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop='BC03lr', speclib = NULL, Z = c(5, 5, 5, 5),
z = 0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, cossplit = c(9e+09, 1.3e+10),
dosplit = FALSE)

SMstarp5(burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10, oldmass = 1e+10,
ancientmass = 1e+10, burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09),
midage = c(1e+09, 5e+09), oldage = c(5e+09, 9e+09), ancientage = c(9e+09, 1.3e+10),
stellpop = "BC03lr", speclib = NULL, Z = c(5, 5, 5, 5, 5), z = 0, H0 = 67.8,
OmegaM = 0.308, OmegaL = 1 - OmegaM, cossplit = c(9e+09, 1.3e+10), dosplit = FALSE)

SMstarfunc(massfunc = function(age, SFR = 1) { ifelse(age < 1e+10, SFR, 0)},
forcemass = FALSE, unimax = 1.38e+10, agescale = 1, stellpop = "BC03lr", speclib = NULL,
Z = 5, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burstmass}{
Numeric scalar; the total stellar mass formed in a recent burst.
}
  \item{youngmass}{
Numeric scalar; the total stellar mass formed in young stars.
}
  \item{midmass}{
Numeric scalar; the total stellar mass formed in mid stars.
}
  \item{oldmass}{
Numeric scalar; the total stellar mass formed in old stars.
}
  \item{ancientmass}{
Numeric scalar; the total stellar mass formed in ancient stars.
}
  \item{burstage}{
Numeric vector; the age limits in years for the burst phase.
}
  \item{massfunc}{
Function; a function which will take \option{age} (in Yrs) as an input and return SFR (in Msun/Yr). Can be passed further arguments via \dots (this is useful if you do not want to hard code values into a function, i.e. for fitting purposes).
}
  \item{youngage}{
Numeric vector; the age limits in years for the young phase.
}
  \item{midage}{
Numeric vector; the age limits in years for the mid phase.
}
  \item{oldage}{
Numeric vector; the age limits in years for the old phase.
}
  \item{ancientage}{
Numeric vector; the age limits in years for the ancient phase.
}
  \item{forcemass}{
Logical or numeric scalar; if FALSE then the stellar mass is directly computed from the \option{massfunc} function, i.e. the output of this must be a true star formation rate in Msun per year. If \option{forcemass} is a number then this specifies the required integrated stellar mass that must be formed over all time, and the output of \option{massfunc} will linearly scaled as required to achieve this.
}
  \item{unimax}{
Numeric scalar; maximum allowed age of any stellar population relative to z=0 (i.e. today). Any star formation implied by \option{massfunc} that is older than this will be set to 0.
}
  \item{agescale}{
Numeric scalar; scaling to apply to the ages (which will passed in years) to \option{massfunc}. E.g. if left at \option{agescale}=1 then \option{massfunc} is expecting the age in years, but if set to 1e-6/1e-9 then \option{massfunc} is expecting the age in Myrs/Gyrs.
}
  \item{stellpop}{
Character scalar; which stellar population library to use, one of \code{\link{BC03lr}} (default), \code{\link{BC03hr}}, \code{\link{EMILES}}.
}
  \item{speclib}{
Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{tau_birth}{
Numeric scalar; dust tau for birth clouds. Associated with sub 100 Myr star formation.
}
  \item{tau_screen}{
Numeric scalar; dust tau for screen.
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters.
}
  \item{Z}{
Integer vector of length 4; specifies which metallicity to use for each star formation phase.
}
  \item{z}{
Numeric scalar; cosmological redshift, where z must be > -1.
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=67.8 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.308).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1 - OmegaM = 0.692).
}
  \item{outtype}{
Character scalar; specifies the type of photometry computed, must be one of 'mag' / 'magAB' (for AB magnitudes, the default), 'Jansky' / 'Jy' (for Jansky fluxes) or 'CGS' (for CGS fluxes). See \code{\link{magABcalc}} for details.
}
  \item{cossplit}{
Numeric vector; the cosmological age at which the ancient star formation event begins and ends (only used if \option{dosplit=TRUE}).
}
  \item{dosplit}{
Boolean scalar; should the \option{cossplit} ancient star formation limits be used?
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Boolean scalar; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.  
}
  \item{\dots}{
Further arguments to be passed to \code{massfunc} function.
}
}
\details{
Optionally the definition of old and ancient can be associated with a fixed time the history of the Universe (e.g. the peak of starformation). This means the exact lookback period adapts as a function of the redshift provided. This often makes more sense then simply fixing the age of the ancient population since it precents the stellar populations becoming older than the Universe.

We use the simple dust models for birth clouds and dust screens as given by Charlot and Fall (2000):

\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}

The defaults should be reasonable in the regime they are used over, where birth cloud dust should only be applied to stellar populations younger than 10 Myrs. Note that BC03 has 70 stellar population spectra which are younger than this, whilst EMILES has none. This means that EMILES cannot realistically capture the birth dust attenuation phase, so for highly star bursting populations BC03 should probably be preferred.
}
\value{
\code{SFHp4} / \code{SFHp5} / \code{SFHfunc} returns a list containing:

\item{flux}{Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux (erg/s/cm^2/Angstrom).}
\item{lum}{Numeric vector; intrinsic luminosity in the standard BC03 units of Lsun/Angstrom.}
\item{out}{Data frame; first column is the filter name, second columns is the central wavelength, third column is the requested output (e.g. AB mag).}
\item{massvec}{Numeric vector; the star formation rate (Msun/year) at each age sampled.}
\item{masstot}{Numeric scalar; total mass formed in solar masses.}
\item{lumtot}{Numeric scalar; total luminosity in solar luminosities.}
\item{M2L}{Numeric scalar; total mass to light ratio.}
\item{call}{The original function call made.}
\item{ages}{Numeric Matrix; age info where rows are the four phases of star formation, columns are the start [1], end [2], duration [3] and mean time [4] of the phases in years.}
\item{masses}{Numeric Matrix; mass info where rows are the four phases of star formation, columns are the total formed in that phase [1] and the total formed over all time [2].}
\item{SFR}{Numeric vector; the star formation rate averaged out over each star formation phase.}
\item{sSFR}{Numeric vector; the speicifc star formation rate averaged out over each star formation phase.}

\code{SMstarp4} / \code{SMstarp5} / \code{SMfunc} is a 5/6/1 element vector containing the stellar mass observed as produced by each phase of star formation (for \code{SMstarp4} and \code{SMstarp5}), and the total stellar mass formed (\code{SMstarp4} / \code{SMstarp5} / \code{SMfunc}, this is the last element listed).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{photom}}, \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}
}
\examples{
SFHdemo_dust=SFHp4()
plot(SFHdemo_dust$out[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)

SMdemo=SMstarp4()
SMdemo

SFHdemo_nodust=SFHp4(tau_birth=0, tau_screen=0)
wave_star=SFHdemo_nodust$flux[,1]

total_atten=sum(c(0,diff(wave_star))*(SFHdemo_nodust$flux[,2]-SFHdemo_dust$flux[,2]))

plot(SFHdemo_nodust$flux, log='xy', ylim=c(1e-20,1e-15), type='l',
xlab=BC03lr$Labels$Wavelab, ylab='Flux (erg/s/cm^2/Ang)')
lines(SFHdemo_dust$flux, col='grey')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), z=0.1, norm=total_atten),
col='brown')
Dale_temp=Dale_interp(type='NormSFR')
lines(Dale_Msol$Wave, Dale_temp$Aspec*total_atten, col='red')

wave=sort(c(SFHdemo_dust$flux[,'wave'], Dale_Msol$Wave))
SpecAtten=approxfun(SFHdemo_dust$flux, rule=2)(wave)+
          approxfun(Dale_Msol$Wave, Dale_temp$Aspec, rule=2)(wave/1.1)
lines(wave, SpecAtten, col='darkgreen')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
