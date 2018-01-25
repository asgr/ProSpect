\name{dustmass}
\alias{dustmass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Dust Mass
}
\description{
Given an input attenuated and unattenuated galaxy spectra and a dust spectrum, calculates the amount of mass in dust.
}
\usage{
dustmass(wave_star, flux_star_nodust, flux_star_dust, wave_dust, lum_dust, z = 0.1,
H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0,
ref = "planck")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave_star}{
Numeric vector; observed stellar wavelength in Angstroms (same number of elements as \option{flux_star_nodust} and \option{flux_star_dust}.
}
  \item{flux_star_nodust}{
Numeric vector; predicted unattenuated stellar flux (erg/s/cm^2/Angstrom).
}
  \item{flux_star_dust}{
Numeric vector; predicted attenuated stellar flux (erg/s/cm^2/Angstrom).
}
  \item{wave_dust}{
Numeric vector; observed dust wavelength in Angstroms (same number of elements as \option{lum_dust}).
}
  \item{lum_dust}{
Numeric vector; observed luminosity in units Lsol / Ang (for 1 Msol of dust), as per \code{\link{blackbody}}, \code{\link{greybody}} and \code{\link{Dale_Msol}}.
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
This calculates the self consistent mass in dust for a given combination of an unattenuated and attenuated spectrum for a given dust template.
}
\value{
Numeric vector; contains the dust mass in solar units (DustMass [1]), dust luminsity in solar units (DustLum [2]) and the dust mass to luminosity ratio M2L [3].
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{blackbody}}, \code{\link{greybody}}, \code{\link{Dale_Msol}}
}
\examples{
#Let's try varying our burst mass and seeing what happen to various outputs:
#The below spans the typical DustMass/StellMass = 0.001 - 0.01 range nicely!

burstmass=0
SFHdemo_dust=SFHp4(burstmass=burstmass)
SFHdemo_nodust=SFHp4(burstmass=burstmass, tau_birth=0, tau_screen=0)
SMdemo=SMstarp4(burstmass=burstmass)
SMdemo
dustout=dustmass(SFHdemo_dust$flux[,1], SFHdemo_nodust$flux[,2], SFHdemo_dust$flux[,2],
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[5]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_nodust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[5]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

burstmass=1e8
SFHdemo_dust=SFHp4(burstmass=burstmass)
SFHdemo_nodust=SFHp4(burstmass=burstmass, tau_birth=0, tau_screen=0)
SMdemo=SMstarp4(burstmass=burstmass)
SMdemo
dustout=dustmass(SFHdemo_dust$flux[,1], SFHdemo_nodust$flux[,2], SFHdemo_dust$flux[,2],
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[5]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_nodust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[5]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

burstmass=1e9
SFHdemo_dust=SFHp4(burstmass=burstmass)
SFHdemo_nodust=SFHp4(burstmass=burstmass, tau_birth=0, tau_screen=0)
SMdemo=SMstarp4(burstmass=burstmass)
SMdemo
dustout=dustmass(SFHdemo_dust$flux[,1], SFHdemo_nodust$flux[,2], SFHdemo_dust$flux[,2],
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[5]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_nodust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[5]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

burstmass=1e10
SFHdemo_dust=SFHp4(burstmass=burstmass)
SFHdemo_nodust=SFHp4(burstmass=burstmass, tau_birth=0, tau_screen=0)
SMdemo=SMstarp4(burstmass=burstmass)
SMdemo
dustout=dustmass(SFHdemo_dust$flux[,1], SFHdemo_nodust$flux[,2], SFHdemo_dust$flux[,2],
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[5]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_nodust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[5]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
