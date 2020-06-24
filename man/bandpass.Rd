\name{bandpass}
\alias{bandpass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Band Pass
}
\description{
Put a spectrum through a filter response, returning either the modified spectrum or the integrated flux.
}
\usage{
bandpass(wave, flux, filter, flux_in = 'freq', flux_out = 'freq',
detect_type = 'photon')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Observed wavelength in Angstroms by default (see \option{wavefac}). For convenience, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{flux}.
}
  \item{flux}{
Numeric vector; flux spectral density in erg/s / cm^2 / Hz or erg/s / cm^2 / Ang units (see \option{flux_in} and \option{flux_out} and conversion in the Example).
}
  \item{filter}{
Two column matrix; column one must be flux in Angstroms, column two must be a filter response. These are the type returned by \code{\link{getfilt}} (see Examples). You can also pass in a function (e.g. the output of \code{\link{approxfun}} run on the filters), which will result in faster processing.
}
  \item{flux_in}{
Character scalar; the input type of flux spectral density used. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
  \item{flux_out}{
Character scalar; the output type of flux spectral density desired. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
  \item{detect_type}{
Character scalar; the type of detector used for the photometric system. Should be either 'photon' (default) or 'energy'. See Details.
}
}
\details{
This is a lower level function that might be useful to use directly. For most use cases a direct conversion using e.g. \code{\link{magABcalc}} is probably more useful. \cr

The type of detector is important to ensure the correct sort of filter operation is made to calculate the average flux at the effective wavelength.

For photon counters we must use (where R is the filter response/transmission):

\deqn{
F(\lambda_{eff}) = \frac{\int R(\lambda) F(\lambda) \lambda d\lambda}{\int R(\lambda) \lambda d\lambda}
}{
F_eff=Int(R_l.F_l.l.dl)/Int(R_l.l.dl)
}

For energy counters we must use:

\deqn{
F(\lambda_{eff}) = \frac{\int R(\lambda) F(\lambda) d\lambda}{\int R(\lambda) d\lambda}
}{
F_eff=Int(R_l.F_l.dl)/Int(R_l.dl)
}

Note the \eqn{F(\lambda) \lambda d\lambda}{F_l.l.dl} parts can be directly replaced by \eqn{F(\nu) \lambda d\nu}{F_n.l.dn}.

}
\value{
Returns the mean weighted flux across the filter.
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
bandpass(fluxBC03[,1], fluxBC03[,2], flux_in='wave', filter=getfilt('r_VST'))

#Alternatively we can output in <F_lambda> and use the pivot wavelength of r_VST
#(6312.451 Ang) to convert to <F_freq>:

bandpass(fluxBC03[,1], fluxBC03[,2], flux_in='wave', flux_out='wave',
filter=getfilt('r_VST')) * (6312.451^2)*1e-10/299792458
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ photometry }% use one of  RShowDoc("KEYWORDS")
