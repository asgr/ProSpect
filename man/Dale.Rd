\name{Dale}
\alias{Dale}
\alias{Dale_Orig}
\alias{Dale_Msol}
\alias{Dale_NormTot}
\alias{Dale_NormAGN}
\alias{Dale_NormSFR}
\alias{Dale_MsolM2L}
\alias{Dale_scale}
\docType{data}
\title{
Dale 2014 Infrared Dust Templates
}
\description{
The full compendium of Dale et al (2014) dust templates scaled for different use cases.

Dale 2014 provides templates that include star formation (SFR) and AGN components. Care needs to be taken when applying them for this reason.
}
\usage{
data("Dale_Orig")
data("Dale_Msol")
data("Dale_NormTot")
data("Dale_NormAGN")
data("Dale_NormSFR")
data("Dale_MsolM2L")
Dale_scale(alpha_SF = 1.5, AGNfrac = 0.5, Dale_in)
}
\format{
The two parameter model marginalises over AGN fraction and the alpha slope of the star forming population. To be consistent with how the data is provided the different AGN fractions are put into different list elements, and the alpha values are stored in different rows. The columns represent fluxes in various units, as per the Dale type used and can be matched against the corresponding length Wave list element.

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
\arguments{
  \item{alpha_SF}{
Numeric scalar; desired interpolated alpha slope of the star forming population.  
}
  \item{AGNfrac}{
Numeric scalar; AGN fraction linearly mixed in energy over the range 5-20 microns (see Dale 2014).
}
  \item{Dale_in}{
Two column matrix; the wavelength and flux of the dust model to be computed. \option{alpha_SF} and \option{AGNfrac} should still be provided correctly, but passing in \option{Dale_in} will often be faster because there will be less data loading required.  
}
}
\details{
A simplified form of the Dale 2014 templates which varies the dust emission by the star formation radiation field (alpha_SF) and AGN intensity (Aspec).

There are 5 different scalings provided:

Dale_Orig: The original W/[H atom] flux scaling, as provided in Dale et al (2014).

Dale_Msol: Modifed flux where units become Lsol / Ang (for 1 Msol of dust). I.e. Dale_Orig$flux.(1.989e30/1.674e-27/3.828e26/0.0073)/Dale_Orig$wave/(1-AGNfrac) [Dale_Msol$Aspec[[i]][j,]=Dale_Orig$Aspec[[i]][j,].(Msol_kg/MH_inkg/Lsol_inW/(M_dust/M_H))/Dale_Orig$wave/Dale_Orig$AGNfrac[i]]. This version can in principle be used to link SFR to FIR directly and self-consistently. M_dust/M_H=0.0073 is taken from Draine (2007) for the MW. The key here is that the scaling is always per unit mass of dust, with AGN flux mixed linearly on top. This means it is easy to fit for dust and AGN and still get back the correct dust mass. See Examples.

Dale_NormTot: Dale_Msol total flux is scaled so it integrates to 1.

Dale_AGNTot: Dale_Msol AGN contributing flux is scaled so it integrates to 1 (so total integral is larger when including the SFR component).

Dale_SFRTot: Dale_Msol SFR contributing flux is scaled so it integrates to 1 (so total integral is larger when including the AGN component).

Dale_MsolM2L: a special output for Dale_Msol. It is a 21 (row) x 64 (column) matrix containing the total dust mass to luminosity ratios for all 21 variants of AGNfrac and 64 variants of alpha_SF respectively. It can be used to infer dust masses for different fitted templates if the known intergral of luminosity is know from the optical dust attenuation. Notice all values are below 1, since dust is much more luminous per unit mass than stars etc (or put a different way, not much dust can attenuate a lot of star light quite effectively, and this energy must be re-radiated).

Something to be aware of is how Dale 2014 does the mixing. It is actually defined as the energy intergral between 5 and 20 microns. This means the "AGNfrac" as computed is not representative of the total energy in AGN compared to attenuated star formation. If you compare our pure SF normalised to 1 versus pure AGN normalised to 1 over this regime (i.e. Dale_NormTot$Aspec[[1]] versus Dale_NormTot$Aspec[[21]]) these actually need to be further scaled by 0.7545396 (for SF) and 0.2454604 (for AGN) to bring them into full alignment. \code{Dale_scale} does this energy scaling for the user, also see examples.
}
\value{
\code{Dale_scale} returns the two component vector bolometric luminosity scaling for dust (element 1, Dustfrac_bol) and AGN (element 2, AGNfrac_bol) components of the FIR emission.
}
\source{
Downloaded the file spectra.SFAGN.tar.gz from http://physics.uwyo.edu/~ddale/research/seds/seds.html.
}
\references{
Dale et al, 2014, ApJ, 784, 11; Draine et al, 2007, ApJ, 663, 866
}
\examples{
### Sanity check of different scaling options:

plot(Dale_Orig$Wave, Dale_Orig$Aspec[[2]][24,], type='l',log='xy')
lines(Dale_Orig$Wave, Dale_Orig$Aspec[[20]][24,],col='red')
lines(Dale_Orig$Wave, Dale_Orig$Aspec[[11]][24,],col='green')

plot(Dale_Msol$Wave, Dale_Msol$Aspec[[2]][24,], type='l',log='xy', ylim=c(1e-12,1))
lines(Dale_Msol$Wave, Dale_Msol$Aspec[[20]][24,],col='red')
lines(Dale_Msol$Wave, Dale_Msol$Aspec[[11]][24,],col='green')

plot(Dale_NormTot$Wave, Dale_NormTot$Aspec[[2]][24,], type='l',log='xy', ylim=c(1e-16,1e-4))
lines(Dale_NormTot$Wave, Dale_NormTot$Aspec[[20]][24,],col='red')
lines(Dale_NormTot$Wave, Dale_NormTot$Aspec[[11]][24,],col='green')

plot(Dale_NormAGN$Wave, Dale_NormAGN$Aspec[[2]][24,], type='l',log='xy', ylim=c(1e-16,1e-4))
lines(Dale_NormAGN$Wave, Dale_NormAGN$Aspec[[20]][24,],col='red')
lines(Dale_NormAGN$Wave, Dale_NormAGN$Aspec[[11]][24,],col='green')

plot(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[2]][24,], type='l',log='xy', ylim=c(1e-16,1e-4))
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[20]][24,],col='red')
lines(Dale_NormSFR$Wave, Dale_NormSFR$Aspec[[11]][24,],col='green')

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

#Check our Dale_scale function for an interpolated AGNfrac and alpha_SF:

alpha_SF=1.53
AGNfrac=0.57

Dale_mix=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, type='NormTot')
Dale_AGN=Dale_interp(alpha_SF=alpha_SF, AGNfrac=1, type='NormTot')
Dale_SFR=Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, type='NormTot')

tempscale=Dale_scale(alpha_SF=alpha_SF, AGNfrac=AGNfrac)

plot(Dale_NormTot$Wave, Dale_mix$Aspec, type='l', log='xy', lwd=4)
lines(Dale_NormTot$Wave, Dale_SFR$Aspec*tempscale[1], col='red')
lines(Dale_NormTot$Wave, Dale_AGN$Aspec*tempscale[2], col='brown')
lines(Dale_NormTot$Wave,
Dale_SFR$Aspec*tempscale[1] + Dale_AGN$Aspec*tempscale[2],
col='green', lty=3, lwd=2)

#The above is how you would reconstruct a generic fit, where you need to scale these
#numbers by the total attenuated flux/luminosity (since above is normalised to 1).

}
\keyword{datasets}
