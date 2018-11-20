\name{CF_general}
\alias{CF_general}
\alias{CF_birth}
\alias{CF_screen}
\alias{CF_birth_atten}
\alias{CF_screen_atten}
\alias{CF}
\alias{dust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Charlot and Fall Simple Dust Model
}
\description{
Applies a simple Charlot and Fall dust model, assuming the published defaults for the birth cloud and dust screen exponential taus.
}
\usage{
CF_birth(wave, tau = 1.0, pow = -0.7, pivot = 5500)
CF_screen(wave, tau = 0.3, pow = -0.7, pivot = 5500)
CF_birth_atten(wave, flux, tau = 1.0, pow = -0.7, pivot = 5500)
CF_screen_atten(wave, flux, tau = 0.3, pow = -0.7, pivot = 5500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; wavelength in Angstroms.
}
  \item{flux}{
Numeric vector; flux in erg/s / cm^2 / Angstrom specral units.
}
  \item{tau}{
Numeric scalar; the exponential tau of the dust, where larger means more reddening. Charlot and Fall suggested defaults are used here.
}
  \item{pow}{
Numeric scalar; power to further raise the attenuation. Charlot and Fall suggest this should be close to -0.7 based on observations.
}
  \item{pivot}{
Numeric scalar; the pivot wavelength in Angstroms. Charlot and Fall suggested default.
}
}
\details{
We use the simple dust models for birth clouds and dust screens as given by Charlot and Fall (2000):

\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}

The defaults should be reasonable in the regime they are used over, where birth cloud dust should only be applied to stellar populations younger than 10 Myrs. Note that BC03 has 70 stellar population spectra which are younger than this, whilst EMILES has none. This means that EMILES cannot realistically capture the birth dust attenuation phase, so for highly star bursting populations BC03 should probably be preferred.
}
\value{
For CF_birth and CF_screen, a numeric vector; the attenuation curve that fluxes should be multiplied by.

For CF_birth_atten and CF_screen_atten, a list where the first element is the attenuated flux [flux], the second is the total energy of the flux attenuated (i.e. this is the amount that should be re-radiated elsewhere) [total_atten], and the third is the ratio of attenuated to unattenuated stellar light (i.e. 1 would mean no flux has been attenuated by dust, and 0 would mean 100\% attenuation).
}
\references{
Charlot & Fall, 2000, ApJ, 539, 718
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
