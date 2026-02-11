\name{greybody}
\alias{blackbody}
\alias{blackbody_norm}
\alias{greybody}
\alias{greybody_norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Blackbody and Greybody Spectra
}
\description{
Produces diferent types of blackbody and greybody spectra.
}
\usage{
blackbody(wave, Temp = 50, k850 = 0.077)
blackbody_norm(wave, Temp = 50, z = 0, norm = 1)
greybody(wave, Temp = 50, beta = 1.5, k850 = 0.077)
greybody_norm(wave, Temp = 50, beta = 1.5, z = 0, norm = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Observed wavelength in Angstroms.
}
  \item{Temp}{
Numeric vector; the absolute temperature of the system in Kelvin (K). For \code{blackbody_norm} and \code{greybody_norm}, if more than one value is provided then multiple black/greybodies are summed together to make a composite greybody spectrum.
}
  \item{z}{
Numeric scalar; cosmological redshift, where z must be > 0.
  }
  \item{beta}{
Numeric vector; represents emissivity, or the spectral emissivity index. This is usually in the range 1-2 (see Hildebrand 1983). For \code{blackbody_norm} and \code{greybody_norm}, this should have the same number of elements as provided for \option{Temp}, but if length=1 then the values are recycled.
}
  \item{k850}{
Numeric scalar; the dust mass normalisation. The default is as used in MAGPHYS (0.077 m^2/kg). A sensible alternative might be 0.0383 m^2/kg (Draine 2003). The latter requires almost twice as much dust to achieve the same luminosity (see Clemens 2013 which gives the general warning "Absolute values for dust masses are, in any case, to be used with caution.").
}
  \item{norm}{
Numeric vector; for \code{blackbody_norm} and \code{greybody_norm}, the integral of the returned black/greybody spectrum over all wavelength. By default it integrates to 1. This should have the same number of elements as provided for \option{Temp}, but if length=1 then the values are recycled.
}
}
\details{
We use a modified version of Eqn 2 from Casey et al (2012). The idea is the normalised output can be scaled to equal any amount of required energy in a simple manner, i.e. if the difference between intrinsic and attenuated flux is known then this can be trivially mapped onto output greybody dust by scaling by this same factor using the \option{norm} argument (since both sets of fluxes are per angstrom).
}
\value{
Numeric vector; observed luminosity in units Lsol / Ang (for 1 Msol of dust) for \code{blackbody} and \code{greybody}, and normalised as specified by \option{norm} for \code{blackbody_norm} and \code{greybody_norm}.
}
\references{
Casey C.M., 2012, MNRAS, 425, 3094 \cr
Clemens M.S. et al., 2013, MNRAS 433, 695 \cr
da Cunha E. et al, 2008, MNRAS 388, 1595 \cr
Draine B.T., 2003, ARA&A, 41, 241 \cr
Hildebrand R.H., 1983, QJRAS, 24, 267
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{dust}}, \code{\link{Dale}}
}
\examples{
#The greybody_norm function can take vector inputs:

plot(greybody_norm(10^seq(5,7,by=0.01), norm=1), type='l', log='x', col='grey',
ylim=c(0, 3e-6))
lines(greybody_norm(10^seq(5,7,by=0.01), Temp=100, norm=0.5), col='grey')
lines(greybody_norm(10^seq(5,7,by=0.01), Temp=c(50,100), norm=c(1,0.5)))

### Check against greybody template:

fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
birthBC03_atten=CF_birth_atten(fluxBC03[,1],fluxBC03[,2])
screenBC03_atten=CF_screen_atten(fluxBC03[,1],fluxBC03[,2])
totalatten=birthBC03_atten$total_atten+screenBC03_atten$total_atten

# Here we show 3 pure greybodies, and then 3 Dale templates normalised by SFR with
# differing amounts of AGN contribution. Note even the pure SFR template is not quite
# greybody (there is an excess in the MIR and very-FIR).

plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
type='l', col='red', xlim=c(1e2,1e7))
lines(fluxBC03[,1], (birthBC03_atten$flux+screenBC03_atten$flux)/2, col='black')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1, norm=totalatten),
col='brown', lty=1)
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=1.5, norm=totalatten),
col='brown', lty=2)
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), beta=2, norm=totalatten),
col='brown', lty=3)
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[1]][24,]*totalatten, col='purple')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[11]][24,]*totalatten, col='orange')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[20]][24,]*totalatten, col='green')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dust }% use one of  RShowDoc("KEYWORDS")
\keyword{ greybody }% use one of  RShowDoc("KEYWORDS")
