\name{SFHp4}
\alias{SFHp4}
\alias{SMstarp4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Four Phase Star Formation History
}
\description{
Computes output galaxy spectra using four distinct phases of constant star formation: recent burst, young, old and ancient. Uses either BC03 or EMILES libraries. Each star formation phase can have a different metallicity associated with it (default solar).
}
\usage{
SFHp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop = "BC03lr", tau_birth = 1, tau_screen = 0.3,
filters = "all", Z = c(5, 5, 5, 5), z = 0.1, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = "planck",
outtype = "mag", cossplit = c(9e+09, 1.3e+10), dosplit = FALSE)

SMstarp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop='BC03lr', Z = c(5, 5, 5, 5), z = 0, H0 = 100,
OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0,
ref = 'planck', cossplit = c(9e+09, 1.3e+10), dosplit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burstmass}{
Numeric scalar; the total stellar mass formed in a recent burst.
}
  \item{youngmass}{
Numeric scalar; the total stellar mass formed in young stars.
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
  \item{youngage}{
Numeric vector; the age limits in years for the young phase.
}
  \item{oldage}{
Numeric vector; the age limits in years for the old phase.
}
  \item{ancientage}{
Numeric vector; the age limits in years for the ancient phase.
}
  \item{stellpop}{
Character scalar; which stellar population library to use, one of \code{\link{BC03lr}} (default), \code{\link{BC03hr}}, \code{\link{EMILES}}.
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
  \item{outtype}{
Character scalar; the type of output to generate. Options are mag, cgs and janksky. See \code{\link{magABcalc}} for details.
}
  \item{cossplit}{
Numeric vector; the cosmological age at which the ancient star formation event begins and ends (only used if \option{dosplit=TRUE}).
}
  \item{dosplit}{
Boolean scalar; should the \option{cossplit} ancient star formation limits be used?
}
}
\details{
Optionally the definition of old and ancient can be associated with a fixed time the history of the Universe (e.g. the peak of starformation). This means the exact lookback period adapts as a function of the redshift provided. This often makes more sense then simply fixing the age of the ancient population since it precents the stellar populations becoming older than the Universe.

We use the simple dust models for birth clouds and dust screens as given by Charlot and Fall (2000):

\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}

The defaults should be reasonable in the regime they are used over, where birth cloud dust should only be applied to stellar populations younger than 10 Myrs. Note that BC03 has 70 stellar population spectra which are younger than this, whilst EMILES has none. This means that EMILES cannot realistically capture the birth dust attenuation phase, so for highly star bursting populations BC03 should probably be preferred.
}
\value{
SFHp4 returns a list containing:

\item{flux}{Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux (erg/s/cm^2/Angstrom).}
\item{lum}{Numeric vector; intrinsic luminosity in the standard BC03 units of Lsun/Angstrom.}
\item{out}{Data frame; first column is the filter name, second columns is the central wavelength, third column is the requested output (e.g. AB mag).}
\item{masstot}{Numeric scalar; total mass formed in solar masses.}
\item{lumtot}{Numeric scalar; total luminosity in solar luminosities.}
\item{M2L}{Numeric scalar; total mass to light ratio.}
\item{call}{The original function call made.}
\item{ages}{Numeric Matrix; age info where rows are the four phases of star formation, columns are the start [1], end [2], duration [3] and mean time [4] of the phases in years.}
\item{masses}{Numeric Matrix; mass info where rows are the four phases of star formation, columns are the total formed in that phase [1] and the total formed over all time [2].}
\item{SFR}{Numeric vector; the star formation rate averaged out over each star formation phase.}
\item{sSFR}{Numeric vector; the speicifc star formation rate averaged out over each star formation phase.}

SMstarp4 is a 5 element vector containing the stellar mass observed as produced by each phase of star formation, and thte total stellar mass formed.
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
total_atten=sum(c(0,diff(BC03lr$Wave))*(SFHdemo_nodust$flux-SFHdemo_dust$flux))

plot(SFHdemo_nodust$flux, log='xy', ylim=c(1e-20,1e-15), type='l',
xlab=BC03lr$Labels$Wavelab, ylab='Flux (erg/s/cm^2/Ang)')
lines(SFHdemo_dust$flux, col='grey')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), z=0.1, norm=total_atten),
col='brown')
Dale_temp=Dale_interp(type='NormSFR')*total_atten
lines(Dale_Msol$Wave, Dale_temp*total_atten, col='red')

wave=sort(c(SFHdemo_dust$flux[,'wave'], Dale_Msol$Wave))
SpecAtten=approxfun(SFHdemo_dust$flux, rule=2)(wave)+
          approxfun(Dale_Msol$Wave,Dale_temp, rule=2)(wave)
lines(wave, SpecAtten, col='darkgreen')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
