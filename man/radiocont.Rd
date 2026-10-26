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
radiocont(wave, flux, z = 0, Te = 10000, ff_frac = 0.1, ff_power = -0.1, 
    sy_power = -0.8, wavesamp = seq(6, 9.4, by = 0.1), flux_in = "freq", 
    flux_out = flux_in, subtractonly = FALSE) 
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
Numeric scalar; effective HII nebular plasma temperature in Kelvin. The free-free flux scales as (Te/1e4)^0.45 at fixed FIR (Marvil et al. 2015, Eq. 4), so this rescales the whole radio SED. Leave as the default (10,000 K), unless you know what you are doing.
}
  \item{ff_frac}{
Numeric scalar; free-free HII nebular plasma thermal radio emission contribution fraction at 1.4 Ghz, i.e. S_ff/S_tot at 1.4 Ghz (Marvil et al. 2015, Eq. 5). Should be between 0.05 and 0.2 typically. By defintion the radio synchrotron emission fraction at 1.4 Ghz will be 1 - \option{ff_frac}.
}
  \item{ff_power}{
Numeric scalar; power-law slope of the free-free nebular plasma thermal radio emission. Leave as the default (-0.1), unless you know what you are doing.
}
  \item{sy_power}{
Numeric scalar; power-law slope of the synchrotron radio emission (should be between -0.8 and -0.6 typically).
}
  \item{wavesamp}{
Numeric vector; desired output log10 wavelength radio emission grid to use in Ang. The default covers the FIR the 1.4 Ghz range, so is usually appropriate for the radio SED extension. If the supplied grid does not reach as long as \code{max(wave)}, it is automatically extended at its native spacing out to \code{max(wave)}. This is necessary because the subtraction always runs over all wavelengths beyond 1e6 Ang, so emission added over a shorter range would leave a hole (i.e. a sharp discontinuity) at the end of the output spectrum.
}
  \item{flux_in}{
Character scalar; the input type of flux spectral density used. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
  \item{flux_out}{
Character scalar; the output type of flux spectral density desired. Should be either 'freq' (default, power per Hz) or 'wave' (power per Ang)
}
  \item{subtractonly}{
Logical; if true then only subtracts the Dale template radio emission, if false also adds on new radio emission. 
}
}
\details{
This radio continuum extension for \code{ProSpect} is mostly based on the work of Marvil et al (2015), and follows their equations directly.

The thermal (free-free) term is their Eq. 4, which gives the free-free flux density in Jy from the FIR flux: S_th(nu) = 1.4e10 * FIR * (Te/1e4)^0.45 * (nu/GHz)^-0.1, for FIR in W/m^2. The code evaluates this at 1.4 Ghz and applies \option{ff_power} as the frequency slope, i.e. S_ff(nu) = 1.4e10 * FIR * (Te/1e4)^0.45 * (nu/1.4Ghz)^ff_power. Note the (Te/1e4)^0.45 dependence is the one appropriate to fixed FIR. This is consistent with their Eq. 2 (SFR proportional to L_th * (Te/1e4)^-0.45), which is written solved for the SFR; inverting Eq. 2 to get L_th at fixed SFR gives the same (Te/1e4)^+0.45. The negative exponent in Eq. 2 as printed is therefore not a sign discrepancy with Eq. 4, and heating the HII gas raises both L_th at fixed SFR and S_th at fixed FIR. The synchrotron term is then fixed by \option{ff_frac} rather than derived independently, via the thermal fraction definition of their Eq. 5, f_th(nu) = S_th(nu)/S(nu).

Note the FIR described above is a proxy for the FIR luminosity rather than an independent measurement of it: it is evaluated from the input spectrum's own flux at 1.4 Ghz, so the radio is being regenerated from the radio already present rather than predicted from the dust emission. This is why the default options reproduce the Dale template radio by construction, but it also means the resulting radio cannot differ from the input radio by more than the change in \option{ff_frac}, \option{Te} and the power-law slopes.

The default options closely reproduce the Dale template radio emission. The advantage of using this interface is the user has fine control over the free-free and synchrotron radio emission mixing (which is fixed to 1:9 at 1.4 Ghz in the Dale templates), and also control over the power-law indices (fixed to -0.1 and -0.8 respectively in the Dale templates). Note the Dale templates are free at exactly \option{ff_frac} = 0.1, so it is that value which reproduces the template radio; other values generate a different, but physically consistent, radio SED.

Condon & Yin (1990) suggest that \option{ff_frac} is seen to typically span the range 0.05 to 0.2, so it should probably be kept within these limits for generative and fitting purposes. Marvil et al (2015) measure a 1.4 Ghz thermal fraction of 9\% +/- 3\% from IRAS far-IR fluxes, which is where the 0.1 default comes from.

Since the free-free term is anchored to the FIR, \option{ff_frac} changes both the free-free to synchrotron mix and the total 1.4 Ghz flux, the latter scaling as 1/\option{ff_frac}; the free-free flux itself is independent of \option{ff_frac}. This is the same statement as Marvil et al (2015) Eq. 8, ^IR f_th,1.4 = 1.7 * 10^(q - 3.53), which relates the 1.4 Ghz thermal fraction to the Helou et al (1985) FIR/radio ratio q. The code's hard-coded FIR calibration is that equation evaluated at f_th = 0.1 (giving q = 2.29, within 0.01 dex of Eq. 8), so the 0.05 to 0.2 range of \option{ff_frac} corresponds to q spanning about 2.00 to 2.60, i.e. about +/-0.3 dex. Care is needed if \option{ff_frac} is fitted jointly with anything else that sets the radio normalisation, since \option{Te} moves the same normalisation (through (Te/1e4)^0.45) and the two are therefore degenerate.

The FIR luminosity is estimated from the flux at 1.4 Ghz (i.e. 2.14e9 Ang rest frame), which must therefore fall within the supplied \option{wave} range for any radio emission to be computed. If it does not, a warning is issued and no radio emission is either subtracted or added (the rest of the spectrum is returned unchanged, though the grid may still be the union of \option{wave} and \option{wavesamp}). Note this requires the input grid to reach log10 wavelength 9.331 or beyond, which is close to the \code{\link{ProSpectSED}} default of \option{waveout} ending at 9.35; shortening \option{waveout} will silently suppress the radio continuum unless this warning is heeded.
}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux in the units of \option{flux_out}.
}
\references{
Condon & Yin, 1990. ApJ, 357, 97 \cr
Helou, Soifer & Rowan-Robinson, 1985. ApJ, 298, L7 \cr
Marvil, Owen & Eilek, 2015, AJ, 149, 32 (arXiv:1408.6296). Eqs 2, 4, 5 and 8 are used here. \cr
Murphy et al, 2011, ApJ, 737, 67
}
\author{
Aaron Robotham
}
\seealso{
\code{\link[ProSpectData:Dale]{Dale}}
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

testSED = ProSpectSED(z=0.05, AGNlum=0, tau_birth=1, tau_screen=0.3,
                    alpha_SF_birth=0.5, alpha_SF_screen=3, speclib=BC03lr, 
                    Dale=Dale_NormTot, AGN=AGN_UnOb_Sparse, filtout=filtout, 
                    Dale_M2L_func=Dale_M2L_func, returnall=TRUE)

# The defaults almost recreate the Dale template emission:

temp = radiocont(testSED$FinalFlux$wave, testSED$FinalFlux$flux, z=0.05)
plot(temp, log='xy', type='l', ylim=c(1e-6,1), xlab='Wave / Ang',
  ylab='Flux Density / Jy', lwd=3)
lines(testSED$FinalFlux$wave, testSED$FinalFlux$flux, col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ radio }% use one of  RShowDoc("KEYWORDS")
