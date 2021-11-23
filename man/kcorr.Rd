\name{kcorr}
\alias{kcorr}
\alias{kcorr_evo}
\alias{Dcorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k and D Corrections
}
\description{
Functions to compute consistently defined k (spectral shift) and D (dust) corrections.
}
\usage{
kcorr(wave, lum, z, DistMod = NULL, filters = "all", H0 = 67.8, OmegaM = 0.308,
  OmegaL = 1 - OmegaM, prospect = NULL)
  
kcorr_evo(wave, lum, z = 10^seq(-1,1,by=0.1), DistMod = NULL, filters = 'all', H0 = 67.8,
  OmegaM = 0.308, OmegaL = 1 - OmegaM, prospect = NULL)

Dcorr(wave, lum_atten, lum_unatten, filters = "all", prospect = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; the spectral wavelength in Angstroms. For convenience in \code{kcorr}, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{lum}.
}
  \item{lum}{
Numeric vector; the spectral luminosity in Lsol / Angstrom (BC03lr / BC03hr and EMILES format). This should be post all astrophysical processes (e.g. \option{FinalLum} output of \code{\link{ProSpectSED}}).
}
  \item{lum_atten}{
Numeric vector; the spectral luminosity in Lsol / Angstrom (BC03lr / BC03hr and EMILES format). This should be the stellar output post dust attenuation (no re-emission, e.g. \option{StarsAtten} output of \code{\link{ProSpectSED}} or \option{lum_atten} output of \code{\link{SFHfunc}}). 
}
  \item{lum_unatten}{
Numeric vector; the spectral luminosity in Lsol / Angstrom (BC03lr / BC03hr and EMILES format). This should be the stellar output pre dust attenuation (e.g. \option{StarsUnAtten} output of \code{\link{ProSpectSED}} or \option{lum_unatten} output of \code{\link{SFHfunc}}). 
}
  \item{z}{
Numeric scalar/vector; cosmological redshift. This is the redshift you wish to make the k-correction at (need not be the redshift of the actual source, but that is most often what people want to know). For \code{kcorr} this should be a scalar, for \code{kcorr_evo} this should be a vector.
}
  \item{DistMod}{
Numeric scalar/vector; distance modulus used where AbsMag_dust = ApMag_dust - DistMod - kcorr, and DistMod = 5log10(LumDist/Mpc) + 25 in units of mag. If provided for speed reasons (because you are repeating the same redshifts and therefore the same \option{DistMod} for a lot of sources) then it must correspond to the distance modulus of \option{z} (if not, wrong results will be returned). See below for Details. For \code{kcorr} this should be a scalar, for \code{kcorr_evo} this should be a vector.
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters. If not a character then this will be passed onto \code{\link{bandpass}}, and may therefore be a list of filter responses, or even better (faster) a list of filter functions.
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.3).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
}
  \item{prospect}{
List; object of class ProSpectSED. This can be passed in for convenience, in which case various inputs above are extracted from the objects structure. Note if you pass this in and specify an argument, the specified argument takes precedent. This way you can easily pass in your \option{prospect} object whilst e.g. modifying the redshift \option{z} (which you might want to do to check for volume limits).
}
}
\details{
Standard convenience functions for computing survey V-maxes etc.

A pure k-correction just adjusts the SED distortion effects of redshifting an object, therefore:

AbsMag_fin = ApMag_fin - DistMod - kcorr

Where in this context AbsMag_fin and ApMag_fin are the SEDs with stellar light both attenuated by dust and re-emitted by dust.

A D-correction takes the above definition of the AbsMag and removes the effect of dust attenuation, taking you all the way to the SED of the pure stellar emission without dust either attenuating or re-emitting:

AbsMag_unatten = AbsMag_atten - Dcorr

Where in this context AbsMag_unatten is the intrinsic SED of the stellar light and AbsMag_atten is the absolute magnitude with stellar light attenuated by dust (no re-emission).

This version of the D-correction can only create more flux in AbsMag_unatten compared to AbsMag_atten, i.e. it ignores the role of dust re-emission which can be significant at longer wavelengths. Note it probably does not make much sense to compute D-corrections beyond the K/W1/W2 bands (depending on the SED and dust prescription) since you are starting to become dominated by net dust re-emission not stellar attenuation and meaning of a D-correction becomes ambiguous.
}
\value{
\code{kcorr} returns a list with:

\item{z}{Numeric scalar; the redshift used.}
\item{DistMod}{Numeric scalar; the distance modulus used.}
\item{AbsMag}{Numeric vector; the per filter rest frame absolute magnitude (dust attenuated, dust re-emitted).}
\item{ApMag}{Numeric vector; the per filter obseved frame apparent magnitude (dust attenuated, dust re-emitted).}
\item{kcorr}{Numeric vector; the implied per filter k-correction (where kcorr = ApMag - DistMod - AbsMag).}

\code{kcorr_evo} returns a list with:

\item{z}{Numeric vector; the redshift used.}
\item{DistMod}{Numeric vector; the distance modulus used.}
\item{AbsMag}{Numeric vector; the per filter rest frame absolute magnitude (dust attenuated, dust re-emitted).}
\item{ApMag}{Numeric matrix; the per filter obseved frame apparent magnitude (dust attenuated, dust re-emitted) as columns; variation with \option{z} as rows.}
\item{kcorr}{Numeric matrix; the implied per filter k-correction (where kcorr = ApMag - DistMod - AbsMag) as columns; variation with \option{z} as rows}

\code{Dcorr} returns a list with:

\item{AbsMag_atten}{Numeric vector; the per filter rest frame absolute magnitude (dust attenuated, no dust re-emitted).}
\item{AbsMag_unatten}{Numeric vector; the per filter rest frame absolute magnitude (dust unattenuated, no dust re-emitted).}
\item{Dcorr}{Numeric vector; the implied per filter D-correction, where: Dcorr = AbsMag_atten - AbsMag_unatten. Should always be positive!}
}
\references{
Interested users should check the book "The Structure and Evolution of Galaxies" by Steven Phillipps (Section 2.7.2).
}
\author{
Aaron Robotham
}
\note{
There can be some ambiguity with k-corrections whether they do or do not fold in the 2.5log10(1 + z) effect of pure bandpass stretching. The confusion here is because imagine a source remained at the same cosmological distance but started receding at high velocity such that its spectrum became purely doppler stretched. Even if the spectrum was perfectly flat with constant flux per band you would still notice dimming because you are numerically spreading the same amount of flux over more wavelength, so the flux density has to drop. So even without having to correct for any change in SED with a function of wavelength, you will still have this band stretching to worry about. The \code{kcorr} function includes this term in it, so even a completely flat spectrum source will have a non-0 \option{kcorr} (remember simple cosmological source dimming is taken care of with the \option{DistMod} term.
}

\seealso{
\code{\link{ProSpectSED}}
}
\examples{
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}

prospect=ProSpectSED(filtout=filtout, z=0.1)

#For k-corrections we can run:

kcorr(prospect$FinalLum, z=0.1, filters=filtout)
#or for convenience
kcorr(prospect=prospect)

#For D-corrections we can run:
Dcorr(prospect$StarsAtten$wave, prospect$StarsAtten$lum, prospect$StarsUnAtten$lum,
  filters=filtout)
#or for convenience
Dcorr(prospect=prospect)
  
#If you want to know a volume limit then you often want to know z vs the apparent mag:

zvec = seq(0.05,1,by=0.05)
ap_evo={}
for(z in zvec){
  #Below we just vary the observed redshift, all other inputs come from prospect
  ap_evo = rbind(ap_evo, kcorr(z=z, prospect=prospect)$ApMag)
}
colnames(ap_evo) = filters

library(magicaxis)
library(celestial)

magplot(zvec, ap_evo[,'r_SDSS'], type='l', xlab='z', ylab='ap_r', xlim=c(0,1))
abline(h=19.8, lty=2) #the nominal GAMA limit
abline(v=0.155, lty=2)

#For GAMA we know the total survey area (293.82 sqdeg) so can compute the Vmax:

cosvol(area=293.82, zmax=0.155, zmin=0.02, H0=67.8, OmegaM=0.308)

#Note we usually apply a lower z limit too because of the large distance errors at low z.

#You can see the impact of not correctly accounting for the k-correction.
#The black dotted line is the k-correction free implied apparent mag, and the red
#shaded region shows the implied k-correction to compute a correct apparent magnitude.

DistMod = cosdistDistMod(zvec, H0=67.8, OmegaM=0.308)
lines(zvec, -19.68021 + DistMod, lty=3) #-19.68021 is the r-band absolute magnitude above
polygon(c(zvec, rev(zvec)), c(ap_evo[,'r_SDSS'], rev(-19.68021 + DistMod)),
  col=hsv(alpha=0.2), border=NA)
}

\concept{ kcorrection }
\concept{ Dcorrection }
