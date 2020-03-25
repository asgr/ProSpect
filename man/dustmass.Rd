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

SFHdemo_dust=SFHfunc(m1=1)
SMdemo=SMstarfunc(m1=1)
SMdemo
dustout=dustmass(SFHdemo_dust$wave_lum, SFHdemo_dust$lum_unatten, SFHdemo_dust$lum_atten,
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[12]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_dust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[12]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

SFHdemo_dust=SFHfunc(m1=10, m2=10)
SMdemo=SMstarfunc(m1=10, m2=10)
SMdemo
dustout=dustmass(SFHdemo_dust$wave_lum, SFHdemo_dust$lum_unatten, SFHdemo_dust$lum_atten,
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[12]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_dust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[12]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

SFHdemo_dust=SFHfunc(m1=100, m2=10)
SMdemo=SMstarfunc(m1=100, m2=10)
SMdemo
dustout=dustmass(SFHdemo_dust$wave_lum, SFHdemo_dust$lum_unatten, SFHdemo_dust$lum_atten,
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[12]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_dust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[12]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))

SFHdemo_dust=SFHfunc(m1=100, m2=100)
SMdemo=SMstarfunc(m1=100, m2=100)
SMdemo
dustout=dustmass(SFHdemo_dust$wave_lum, SFHdemo_dust$lum_unatten, SFHdemo_dust$lum_atten,
Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,])
dustout
dustout[1]/SMdemo[12]
dustflux=Lum2Flux(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][24,]*dustout[1])
plot(SFHdemo_dust$flux, type='l', log='xy', xlim=c(1e2,1e7), ylim=c(1e-20,1e-12),
xlab='Wavelength / Ang', ylab='Flux (erg/s/cm^2/Ang)', col='darkgrey')
lines(dustflux, col='brown')
lines(addspec(dustflux[,1], dustflux[,2], SFHdemo_dust$flux[,1], SFHdemo_dust$flux[,2]))
legend('topleft', legend=paste('SM =',round(log10(SMdemo[12]),2)))
legend('topright', legend=paste('DM =',round(log10(dustout[1]),2)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ dust }% use one of  RShowDoc("KEYWORDS")
