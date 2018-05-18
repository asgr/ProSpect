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
dustmass(wave_star, lum_star_nodust, lum_star_dust, wave_dust, lum_dust)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave_star}{
Numeric vector; observed stellar wavelength in Angstroms (same number of elements as \option{flux_star_nodust} and \option{flux_star_dust}.
}
  \item{lum_star_nodust}{
Numeric vector; predicted unattenuated stellar luminosity in units Lsol / Ang (for 1 Msol).
}
  \item{lum_star_dust}{
Numeric vector; predicted attenuated stellar luminosity in units Lsol / Ang (for 1 Msol).
}
  \item{wave_dust}{
Numeric vector; observed dust wavelength in Angstroms (same number of elements as \option{lum_dust}).
}
  \item{lum_dust}{
Numeric vector; observed luminosity in units Lsol / Ang (for 1 Msol of dust), as per \code{\link{blackbody}}, \code{\link{greybody}} and \code{\link{Dale_Msol}}.
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
dustout=dustmass(BC03lr$Wave, SFHdemo_nodust$lum, SFHdemo_dust$lum,
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
dustout=dustmass(BC03lr$Wave, SFHdemo_nodust$lum, SFHdemo_dust$lum,
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
dustout=dustmass(BC03lr$Wave, SFHdemo_nodust$lum, SFHdemo_dust$lum,
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
dustout=dustmass(BC03lr$Wave, SFHdemo_nodust$lum, SFHdemo_dust$lum,
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
