\name{Dale}
\alias{Dale}
\alias{Dale_Orig}
\alias{Dale_Msol}
\alias{Dale_NormTot}
\alias{Dale_NormAGN}
\alias{Dale_NormSFR}
\docType{data}
\title{
Dale 2014 Infrared Dust Templates
}
\description{
The full compendium of Dale et al (2014) dust templates scaled for different use cases.

Dale 2014 provides templates that include star formation (SFR) and AGN components. Care needs to be taken when applying them for this reason.
}
\usage{data("Dale_Orig")}
\format{
The two parameter model marginalises over AGN fraction and the alpha slope of the star forming population. To be consistent with how the data is provided the different AGN fractions are put into different list elements, and the alpha values are stored in different rows. The columns represent fluxes in

The wavelength (Wave) is in angstroms.

The flux output is scaled differently for the 5 different scaling schemes (see Details).

An example flux format (using Dale_NormTot) is:

List of 4\cr
 $ AGNfrac : num [1:21] 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...\cr
 $ alpha_SF: num [1:64] 0.0625 0.125 0.1875 0.25 0.3125 ...\cr
 $ Wave    : num [1:1496] 3600 4500 5800 7500 10000 12800 16500 21200 27000 32000 ...\cr
 $ Aspec   :List of 21\cr
  ..$ : num [1:64, 1:1496] 2.23e-09 2.25e-09 2.27e-09 2.30e-09 2.33e-09 ...\cr
  ..$ : num [1:64, 1:1496] 1.44e-06 1.43e-06 1.43e-06 1.42e-06 1.42e-06 ...\cr
  ..$ : num [1:64, 1:1496] 2.98e-06 2.97e-06 2.96e-06 2.95e-06 2.95e-06 ...\cr
  ..$ : num [1:64, 1:1496] 4.66e-06 4.64e-06 4.63e-06 4.61e-06 4.60e-06 ...\cr
  ..$ : num [1:64, 1:1496] 6.48e-06 6.45e-06 6.44e-06 6.42e-06 6.40e-06 ...\cr
  ..$ : num [1:64, 1:1496] 8.46e-06 8.43e-06 8.41e-06 8.39e-06 8.36e-06 ...\cr
  ..$ : num [1:64, 1:1496] 1.06e-05 1.06e-05 1.06e-05 1.05e-05 1.05e-05 ...\cr
  ..$ : num [1:64, 1:1496] 1.30e-05 1.30e-05 1.29e-05 1.29e-05 1.29e-05 ...\cr
  ..$ : num [1:64, 1:1496] 1.57e-05 1.56e-05 1.56e-05 1.55e-05 1.55e-05 ...\cr
  ..$ : num [1:64, 1:1496] 1.86e-05 1.85e-05 1.85e-05 1.84e-05 1.84e-05 ...\cr
  ..$ : num [1:64, 1:1496] 2.19e-05 2.18e-05 2.17e-05 2.17e-05 2.16e-05 ...\cr
  ..$ : num [1:64, 1:1496] 2.55e-05 2.54e-05 2.54e-05 2.53e-05 2.53e-05 ...\cr
  ..$ : num [1:64, 1:1496] 2.97e-05 2.96e-05 2.95e-05 2.95e-05 2.94e-05 ...\cr
  ..$ : num [1:64, 1:1496] 3.44e-05 3.43e-05 3.42e-05 3.42e-05 3.41e-05 ...\cr
  ..$ : num [1:64, 1:1496] 3.99e-05 3.98e-05 3.97e-05 3.96e-05 3.96e-05 ...\cr
  ..$ : num [1:64, 1:1496] 4.62e-05 4.61e-05 4.60e-05 4.60e-05 4.59e-05 ...\cr
  ..$ : num [1:64, 1:1496] 5.37e-05 5.36e-05 5.35e-05 5.35e-05 5.34e-05 ...\cr
  ..$ : num [1:64, 1:1496] 6.27e-05 6.26e-05 6.25e-05 6.24e-05 6.24e-05 ...\cr
  ..$ : num [1:64, 1:1496] 7.36e-05 7.35e-05 7.34e-05 7.34e-05 7.33e-05 ...\cr
  ..$ : num [1:64, 1:1496] 8.72e-05 8.71e-05 8.71e-05 8.71e-05 8.70e-05 ...\cr
  ..$ : num [1:64, 1:1496] 0.000105 0.000105 0.000105 0.000105 0.000105 ...\cr
}
\details{
A simplified form of the Dale 2014 templates which varies the dust emission by the star formation radiation field (alpha_SF) and AGN intensity (Aspec).

There are 5 different scalings provided:

Dale_Orig: The original W/[H atom] flux scaling, as provided in Dale et al (2014).

Dale_Msol: Modifed flux where units become Lsol / Ang (for 1 Msol of dust). I.e. Dale_Orig$flux.(1.989e30/1.674e-27/3.828e26/0.0073)/Dale_Orig$wave [Dale_Orig$flux.(Msol_kg/MH_inkg/Lsol_inW/(M_dust/M_H))/Dale_Orig$wave]. This version can in principle be used to link SFR to FIR directly and self-consistently. M_dust/M_H=0.0073 is taken from Draine (2007) for the MW. See Examples.

Dale_NormTot: Dale_Msol total flux is scaled so it integrates to 1.

Dale_AGNTot: Dale_Msol AGN contributing flux is scaled so it integrates to 1 (so total integral is larger when including the SFR component).

Dale_SFRTot: Dale_Msol SFR contributing flux is scaled so it integrates to 1 (so total integral is larger when including the AGN component).

}
\source{
Downloaded the file spectra.SFAGN.tar.gz from http://physics.uwyo.edu/~ddale/research/seds/seds.html.
}
\references{
Dale et al, 2014, ApJ, 784, 11; Draine et al, 2007, ApJ, 663, 866
}
\examples{
fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)
birthBC03_atten=CF_birth_atten(fluxBC03[,1],fluxBC03[,2])
screenBC03_atten=CF_screen_atten(fluxBC03[,1],fluxBC03[,2])
totalatten=birthBC03_atten$total_atten+screenBC03_atten$total_atten

# Here we show a pure greybody, and then 3 Dale templates normalised by SFR with
# differing amounts of AGN contribution. Note even the pure SFR template is not quite
# greybody (there is an excess in the MIR and very-FIR).

plot(fluxBC03, log='xy', xlab=BC03lr$Labels$Wavelab, ylab='Flux / erg/s/cm^2/Angstrom',
type='l', col='red', xlim=c(1e2,1e7))
lines(fluxBC03[,1], birthBC03_atten$flux, col='black')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01),
norm=birthBC03_atten$total_atten), col='brown')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[1]][24,]*totalatten, col='purple')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[11]][24,]*totalatten, col='orange')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[20]][24,]*totalatten, col='green')
}
\keyword{datasets}
