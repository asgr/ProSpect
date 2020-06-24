\name{radiocont}
\alias{radiocont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Radio Emission Extension
}
\description{
Can take an output \code{ProSpect} spectrum and extend it into the radio regime using a mixture of free-free (thermal) and synchrotron (non-thermal) radio emission. Default values are sensible for typical galaxies, and reproduce the Dale template emission.
}
\usage{
radiocont(wave, flux, z = 0, Te = 10000, ff_frac = 0.1, ff_power = -0.1, sy_power = -0.8,
wavesamp = seq(7, 9.35, by = 0.01), flux_in = "freq", flux_out = flux_in)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; observed wavelength in Angstroms. For convenience, if a two column matrix or data frame is supplied then the first column is taken to be wave and the second is taken to be \option{flux}.
}
  \item{flux}{
Numeric vector; flux spectral density in per Hz or per Ang units (see \option{flux_in} and \option{flux_out}).
}
  \item{z}{
Numeric scalar; the redshift of the input spectrum (usually only relevant for observed frame flux spectra, not rest frame luminosity spectra).
}
  \item{Te}{
Numeric scalar; effective HII nebular plasma temperature in Kelvin. Leave as the default (10,000 K), unless you know what you are doing.
}
  \item{ff_frac}{
Numeric scalar; free-free HII nebular plasma thermal radio emission contribution fraction at 1.4 Ghz (should be between 0.05 and 0.2 typically). By defintion the radio synchrotron emission fraction at 1.4 Ghz will be 1 - \option{ff_frac}.
}
  \item{ff_power}{
Numeric scalar; power-law slope of the free-free nebular plasma thermal radio emission. Leave as the default (-0.1), unless you know what you are doing.
}
  \item{sy_power}{
Numeric scalar; power-law slope of the synchrotron radio emission (should be between -0.8 and -0.6 typically).
}
  \item{wavesamp}{
Numeric vector; desired output log10 wavelength radio emission grid to use in Ang. The default covers the FIR the 1.4 Ghz range, so is usually appropriate for the radio SED extension.
}
  \item{flux_in}{
Character scalar; the input type of flux spectral density used. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
  \item{flux_out}{
Character scalar; the output type of flux spectral density desired. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
}
\details{
This radio continuum extension for \code{ProSpect} is mostly based on the work of Marvil et al (2015).

The default options closely reproduce the Dale template radio emission. The advantage of using this interface is the user has fine control over the free-free and synchrotron radio emission mixing (which is fixed to 1:9 at 1.4 Ghz in the Dale templates), and also control over the power-law indices (fixed to -0.1 and -0.8 respectively in the Dale templates).

Condon & Yin suggest that \option{ff_frac} is seen to typically span the range 0.05 to 0.2, so it should probably be kept within these limits for generative and fitting purposes.
}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux in the units of \option{flux_out}.
}
\references{
Condon & Yin, 1990. ApJ, 357, 97 \cr
Marvil et al, 2015, AJ, 149, 32
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Dale}}
}
\examples{
data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS',
  'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE',
  'W3_WISE', 'W4_WISE', 'P100_Herschel', 'P160_Herschel', 'S250_Herschel' ,
  'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}

testSED=ProSpectSED(z=0.05, AGNlum=0, tau_birth=1, tau_screen=0.3,
                    alpha_SF_birth=0.5, alpha_SF_screen=3, speclib=BC03lr, 
                    Dale=Dale_NormTot, AGN=AGN_UnOb_Sparse, filtout=filtout, 
                    Dale_M2L_func=Dale_M2L_func, returnall=TRUE)

# The defaults recreate the Dale template emission:

temp = radiocont(testSED$FinalFlux$wave, testSED$FinalFlux$flux)
plot(temp, log='xy', type='l', ylim=c(1e-6,1), xlab='Wave / Ang',
  ylab='Flux Density / Jy', lwd=3)
lines(testSED$FinalFlux$wave, testSED$FinalFlux$flux, col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ radio }% use one of  RShowDoc("KEYWORDS")
