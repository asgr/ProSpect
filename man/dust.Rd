\name{Attenuate}
\alias{CF}
\alias{CF_atten}
\alias{CF_birth}
\alias{CF_birth_atten}
\alias{CF_screen}
\alias{CF_screen_atten}
\alias{atten_emit}
\alias{dust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Charlot and Fall Simple Dust Model
}
\description{
Applies a simple Charlot and Fall dust model, assuming the published defaults for the birth cloud and dust screen exponential taus.
}
\usage{
CF(wave, tau = 0.3, pow = -0.7, pivot = 5500)
CF_atten(wave, flux, tau = 0.3, pow = -0.7, pivot = 5500)
CF_birth(wave, tau = 1.0, pow = -0.7, pivot = 5500)
CF_birth_atten(wave, flux, tau = 1.0, pow = -0.7, pivot = 5500)
CF_screen(wave, tau = 0.3, pow = -0.7, pivot = 5500, Eb = 0, L0 = 2175.8,
  LFWHM = 470)
CF_screen_atten(wave, flux, tau = 0.3, pow = -0.7, pivot = 5500, Eb = 0,
  L0 = 2175.8, LFWHM = 470)
atten_emit(wave, flux, tau = 0.3, pow = -0.7, alpha_SF = 1.5, Dale = NULL,
Dale_M2L_func = NULL, waveout = NULL, Eb = 0, L0 = 2175.8, LFWHM = 470) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; wavelength in Angstroms.
}
  \item{flux}{
Numeric vector; flux in any units as long as they are per unit wavelength (e.g. Lsol / Ang).
}
  \item{tau}{
Numeric scalar; the exponential tau of the dust, where larger means more reddening. Charlot and Fall suggested defaults are used here (although tau=0.5 for the screen is also reasonable).
}
  \item{pow}{
Numeric scalar; power to further raise the attenuation. Charlot and Fall suggest this should be close to -0.7 based on observations (da Cunha et al (2008) use a steeper -1.3 for the birth cloud).
}
  \item{pivot}{
Numeric scalar; the pivot wavelength in Angstroms. Charlot and Fall suggested default.
}
  \item{Eb}{
Numeric scalar; Additional attenuation at the 2175.8 Ang dust bump. The default of 0 means this has no effect, MW is near 3.3, and typical galaxies are around \option{Eb}=1.
}
  \item{L0}{
Numeric scalar; location of the 2175.8 Ang dust bump in Angstroms (probably do not adjust this). 
}
  \item{LFWHM}{
Numeric scalar; width of the 2175.8 Ang dust bump in Angstroms (probably do not adjust this). 
}
  \item{alpha_SF}{
Numeric scalar; desired interpolated alpha slope of the star forming population. Lower values mean hotter dust.  
}
  \item{Dale}{
Pass in the Dale dust library directly (should be one of \option{Dale_NormTot}, \option{Dale_NormSFR}).
}
  \item{Dale_M2L_func}{
Function; if dust masses are desired then the user should supply a function to convert \option{alpha_SF} to mass-to-light, e.g. \code{\link{Dale_M2L_func}}
}
  \item{waveout}{
Numeric vector; desired output wavelength grid to use for adding together spectra. See \code{\link{addspec}}.
}
}
\details{
We use the simple dust models for birth clouds and dust screens as given by Charlot and Fall (2000):

\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}

The defaults should be reasonable in the regime they are used over, where birth cloud dust should only be applied to stellar populations younger than 10 Myrs. Note that BC03 has 70 stellar population spectra which are younger than this, whilst EMILES has none. This means that EMILES cannot realistically capture the birth dust attenuation phase, so for highly star bursting populations BC03 should probably be preferred.
}
\value{
For \code{CF}, \code{CF_birth} and \code{CF_screen}, a numeric vector; the attenuation curve that fluxes should be multiplied by.

For \code{CF_atten}, \code{CF_birth_atten} and \code{CF_screen_atten}, a list where the first element is the attenuated flux [flux], the second is the total energy of the flux attenuated (i.e. this is the amount that should be re-radiated elsewhere) [total_atten], and the third is the ratio of attenuated to unattenuated stellar light (i.e. 1 would mean no flux has been attenuated by dust, and 0 would mean 100\% attenuation).

For \code{atten_emit} the input spectrum is attenuated by the specified Charlot and Fall mode, and then re-emitted with the specified Dale dust template. Returns a list containing:

\describe{
  \item{final}{Two column data.frame; Wavelength in Ang, attenuated and re-emitted flux in user units.}
  \item{unatten}{Two column data.frame; Wavelength in Ang, un-attenuated flux in user units.}
  \item{atten}{Two column data.frame; Wavelength in Ang, attenuated flux in user units.}
  \item{emit}{Two column data.frame; Wavelength in Ang, re-emitted flux in user units (i.e. mostly FIR).}
  \item{totalatten}{Integrated total attenuation user units.}
  \item{dustmass}{If \option{Dale_M2L_func} is provided then dustmass in solar masses (Msol).}
}
}
\references{
Charlot & Fall, 2000, ApJ, 539, 718 \cr
da Cunha et al, 2008, MNRAS, 388, 1595
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFH}}, \code{\link{greybody}}
}
\examples{
fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
birthBC03=CF_birth(fluxBC03[,1])*fluxBC03[,2]
screenBC03=CF_screen(fluxBC03[,1])*fluxBC03[,2]

plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux (erg/s/cm^2/Ang)',
type='l', col='red')
lines(fluxBC03[,1], birthBC03, col='black')
lines(fluxBC03[,1], screenBC03, col='grey')
legend('topright', legend=c('Intrinsic', 'Birth Cloud Dust', 'Screen Dust'),
col=c('red', 'black', 'grey'), lty=1)

#We can use the total energy attenuated to re-emit the radiation with a greybody spectrum.

birthBC03_atten=CF_birth_atten(fluxBC03[,1],fluxBC03[,2])

plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
type='l', col='red', xlim=c(1e2,1e7))
lines(fluxBC03[,1], birthBC03, col='black')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), z=0,
norm=birthBC03_atten$total_atten), col='brown')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dust }% use one of  RShowDoc("KEYWORDS")
