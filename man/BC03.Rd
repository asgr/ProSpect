\name{BC03}
\alias{BC03}
\alias{BC03lr}
\alias{BC03hr}
\docType{data}
\title{
GALAXEV Bruzual & Charlot 2003 Models (GALEXev)
}
\description{
This is a reformatted collection of all the key informtion for the low/high resolution BC03 models (technically named GALEXev). This gives enough information to generate generic spectra for a target metallicity and star formation history.
}
\usage{
data("BC03lr")
data("BC03hr")
}
\format{
The highest level list structure contains the following elements:

\describe{
  \item{Z}{Numeric vector of length 6 describing the availability and location of different metallicity stellar populations.}
  \item{Age}{Numeric vector of length 221 giving the mid point of the age in years for the stellar populations.}
  \item{AgeBins}{Numeric vector of length 222 giving the age bin extremes in years of the stellar populations.}
  \item{AgeWeights}{Numeric vector of length 221 giving the total time in years for the stellar population bins.}
  \item{Wave}{Numeric vector of length 1221 giving the spectral wavelength in Angstroms.}
  \item{Labels}{Character vector of length 5 giving handy axis labels.}
  \item{Zpsec}{List containing 6 (for each metallicity) 221x1221 matrices that give the full spectral information.}
  \item{Zevo}{List containing 6 (for each metallicity) 221x5 data frames that contain the evolution information for the stellar material (stellar mass in stars [1] / gas [2] / total [3]; the star formation rate [4] and the stellar mass in remnants [5])}
}

The full low resolution format is:\cr
\cr
List of 8\cr
 $ Z         : num [1:6] 1e-04 4e-04 4e-03 8e-03 2e-02 5e-02\cr
 $ Age       : num [1:221] 0 125893 141254 158489 177828 ...\cr
 $ AgeBins   : num [1:222] 0 62946 133574 149872 168158 ...\cr
 $ AgeWeights: num [1:221] 62946 70627 16298 18287 20519 ...\cr
 $ Wave      : num [1:1221] 91 94 96 98 100 102 104 106 108 110 ...\cr
 $ Labels    :List of 5\cr
  ..$ Zlab      : chr "Metallicity"\cr
  ..$ Agelab    : chr "Time since ZAM / Yrs"\cr
  ..$ Wavelab   : chr "Wavelength / Ang"\cr
  ..$ Lumlab   : chr "Lsol / Ang (for 1 Msol SF)"\cr
  ..$ LumAgelab: chr "Lsol / Ang (for 1 Msol/Yr SFR)"\cr
 $ Zspec     :List of 6\cr
  ..$ : num [1:221, 1:1221] 9.09e-08 9.09e-08 9.09e-08 9.09e-08 9.09e-08 ...\cr
  ..$ : num [1:221, 1:1221] 2.11e-08 2.11e-08 2.11e-08 2.11e-08 2.11e-08 ...\cr
  ..$ : num [1:221, 1:1221] 1.11e-10 1.11e-10 1.11e-10 1.11e-10 1.11e-10 ...\cr
  ..$ : num [1:221, 1:1221] 6.34e-11 6.34e-11 6.34e-11 6.34e-11 6.34e-11 ...\cr
  ..$ : num [1:221, 1:1221] 1.04e-14 1.04e-14 1.04e-14 1.04e-14 1.04e-14 ...\cr
  ..$ : num [1:221, 1:1221] 2.65e-18 2.65e-18 2.69e-18 2.68e-18 2.67e-18 ...\cr
 $ Zevo      :List of 6\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  
The full high resolution format is:\cr
\cr
List of 8\cr
 $ Z         : num [1:6] 1e-04 4e-04 4e-03 8e-03 2e-02 5e-02\cr
 $ Age       : num [1:221] 0 125893 141254 158489 177828 ...\cr
 $ AgeBins   : num [1:222] 0 62946 133574 149872 168158 ...\cr
 $ AgeWeights: num [1:221] 62946 70627 16298 18287 20519 ...\cr
 $ Wave      : num [1:6900] 91 94 96 98 100 102 104 106 108 110 ...\cr
 $ Labels    :List of 5\cr
  ..$ Zlab      : chr "Metallicity"\cr
  ..$ Agelab    : chr "Time since ZAM / Yrs"\cr
  ..$ Wavelab   : chr "Wavelength / Ang"\cr
  ..$ Lumlab   : chr "Lsol / Ang (for 1 Msol SF)"\cr
  ..$ LumAgelab: chr "Lsol / Ang (for 1 Msol/Yr SFR)"\cr
 $ Zspec     :List of 6\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 9.08e-08 9.08e-08 9.08e-08 9.08e-08 ...\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 2.11e-08 2.11e-08 2.11e-08 2.11e-08 ...\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 1.11e-10 1.11e-10 1.11e-10 1.11e-10 ...\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 6.34e-11 6.34e-11 6.34e-11 6.34e-11 ...\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 1.04e-14 1.04e-14 1.04e-14 1.04e-14 ...\cr
  ..$ : num [1:221, 1:6900] 6.90e+03 2.65e-18 2.65e-18 2.69e-18 2.68e-18 ...\cr
 $ Zevo      :List of 6\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  ..$ :'data.frame':	221 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SMgas : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMtot : num [1:221] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:221] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:221] 0 0 0 0 0 0 0 0 0 0 ...\cr
}
\details{
There are 6 metallicities (as given in \option{Z}), 221 age steps per population and 1221/6990 spectral fluxes (low/high resolution). Everything has been computed for a Chabrier IMF and the Padova 1994 isochrones (as recommended in BC03).
}
\source{
http://www.bruzual.org/bc03
}
\references{
BC03/GALEXev: Bruzual and Charlot, 2003, MNRAS, 344, 1000 \cr
Padova 1994 isochrones: Bertelli et al., 1994, A\&AS, 106, 275 \cr
Chabrier IMF: Chabrier, 2003, PASP, 115, 763
}
\examples{
#Low res spectrum of Sun today (Z=0.02, Age=4.75 Gyrs) [nearest to 4.6 Gyrs]:
plot(BC03lr$Wave, BC03lr$Zspec[[5]][160,], log='xy', type='l', xlim=c(2e3,4e4),
ylim=c(1e-8,1e-4), xlab=BC03lr$Labels$Wavelab, ylab=BC03lr$Labels$Lumlab)

#High res spectrum of Sun today (Z=0.02, Age=4.75 Gyrs) [nearest to 4.6 Gyrs]:
plot(BC03hr$Wave, BC03hr$Zspec[[5]][160,], log='xy', type='l', xlim=c(2e3,4e4),
ylim=c(1e-8,1e-4), xlab=BC03lr$Labels$Wavelab, ylab=BC03lr$Labels$Lumlab)
}
\keyword{datasets}
