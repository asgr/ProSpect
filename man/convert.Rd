\name{convert_flux}
\alias{convert_wave2freq}
\alias{convert_freq2wave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Flux From per Wavelength to per Frequency
}
\description{
Conversion routines that by default transform back and forth between spectral flux specified in erg/s / cm^2 / Angstrom units and erg/s / cm^2 / Hz units. The former is the more typical way to plot spectral flux (e.g. galaxy spectra etc), whilst the latter is what is actually needed in order to compute AB magnitudes, CGS flux and Janskys.
}
\usage{
convert_wave2freq(flux_wave, wave, wavefac = 1e-10, freqfac = 1)
convert_freq2wave(flux_freq, wave, wavefac = 1e-10, freqfac = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flux_wave}{
Numeric vector; flux in erg/s / cm^2 / Wavelength spectral flux density units.
}
\item{flux_freq}{
Numeric vector; flux in erg/s / cm^2 / Frequency spectral flux density units.
}
  \item{wave}{
Numeric vector; wavelength in units of \option{wavefac}.
}
  \item{wavefac}{
Numeric scalar; specified unit type converting 1m to chosen unit. Default of 1e-10 means Angstrom, whilst 1e-9 would be nm etc.
}
  \item{freqfac}{
Numeric scalar; specified unit type converting 1Hz to chosen unit. Default of 1 means Hz.
}
}
\value{
Numeric vector; the converted spectral flux.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{magABcalc}}, \code{\link{getfilt}}, \code{\link{SFHfunc}}
}
\examples{
fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)

#The below should give the same answers (6.515128e-28):

CGScalc(fluxBC03)
bandpass(fluxBC03[,1], (1e-10*fluxBC03[,2]*fluxBC03[,1]^2)/299792458,
filter=getfilt('r_VST'))
bandpass(fluxBC03[,1], convert_wave2freq(fluxBC03[,2], fluxBC03[,1]),
filter=getfilt('r_VST'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ flux }% use one of  RShowDoc("KEYWORDS")
