\name{BPASS}
\alias{BPASS}
\docType{data}
\title{
BPASS v2.2.1 Binaries Stellar Library
}
\description{
This is a reformatted collection of all the key informtion for the BPASS v2.2.1 models (binary star, Chabrier IMF up to 300 Msol variant). This gives enough information to generate generic spectra for a target metallicity and star formation history.
}
\usage{data("BPASS")}
\format{
The highest level list structure contains the following elements:

\describe{
  \item{Z}{Numeric vector of length 13 describing the availability and location of different metallicity stellar populations.}
  \item{Age}{Numeric vector of length 51 giving the mid point of the age in years for the stellar populations.}
  \item{AgeBins}{Numeric vector of length 52 giving the age bin extremes in years of the stellar populations.}
  \item{AgeWeights}{Numeric vector of length 51 giving the total time in years for the stellar population bins.}
  \item{Wave}{Numeric vector of length 1221 giving the spectral wavelength in Angstroms.}
  \item{Labels}{Character vector of length 5 giving handy axis labels.}
  \item{Zpsec}{List containing 13 (for each metallicity) 51x1221 matrices that give the full spectral information.}
  \item{Zevo}{List containing 13 (for each metallicity) 51x5 data frames that contain the evolution information for the stellar material (stellar mass in stars [1] / gas [2] / total [3]; the star formation rate [4] and the stellar mass in remnants [5])}
}

The full low resolution format is:\cr
\cr
List of 8 \cr
 $ Z         : num [1:13] 0.00001 0.0001 0.001 0.002 0.003 0.004 0.006 0.008 0.01 0.014 ... \cr
 $ Age       : num [1:51] 1000000 1258925 1584893 1995262 2511886 ... \cr
 $ AgeBins   : num [1:52] 1000000 1122018 1412538 1778279 2238721 ... \cr
 $ AgeWeights: num [1:51] 1122018 290519 365742 460442 579662 ... \cr
 $ Wave      : num [1:1221] 91 94 96 98 100 102 104 106 108 110 ... \cr
 $ Labels    :List of 5 \cr
  ..$ Zlab     : chr "Metallicity" \cr
  ..$ Agelab   : chr "Time since ZAM / Yrs" \cr
  ..$ Wavelab  : chr "Wavelength / Ang" \cr
  ..$ Lumlab   : chr "Lsun / Ang (for 1 Msun SF)" \cr
  ..$ LumAgelab: chr "Lsun / Ang (for 1 Msun/Yr SFR)" \cr
 $ Zspec     :List of 13 \cr
  ..$ : num [1:51, 1:1221] 1.24e-11 1.40e-11 2.27e-11 1.02e-04 4.07e-03 ... \cr
  ..$ : num [1:51, 1:1221] 8.40e-11 9.50e-11 2.03e-10 2.50e-04 1.88e-03 ... \cr
  ..$ : num [1:51, 1:1221] 2.30e-10 3.40e-10 2.91e-10 2.29e-04 1.55e-03 ... \cr
  ..$ : num [1:51, 1:1221] 2.88e-10 2.21e-10 1.49e-10 4.02e-05 3.62e-04 ... \cr
  ..$ : num [1:51, 1:1221] 2.77e-10 1.46e-10 9.46e-11 1.01e-05 2.24e-05 ... \cr
  ..$ : num [1:51, 1:1221] 2.18e-10 9.05e-11 6.24e-11 2.24e-09 7.06e-09 ... \cr
  ..$ : num [1:51, 1:1221] 1.16e-10 2.76e-11 1.80e-11 2.04e-10 1.78e-10 ... \cr
  ..$ : num [1:51, 1:1221] 5.63e-11 1.02e-11 7.65e-12 4.99e-12 2.37e-11 ... \cr
  ..$ : num [1:51, 1:1221] 1.21e-11 2.18e-12 1.47e-12 1.16e-12 9.30e-13 ... \cr
  ..$ : num [1:51, 1:1221] 2.44e-12 1.83e-13 1.69e-13 1.56e-13 1.50e-13 ... \cr
  ..$ : num [1:51, 1:1221] 1.87e-13 1.33e-14 1.29e-14 1.09e-14 1.13e-14 ... \cr
  ..$ : num [1:51, 1:1221] 8.40e-16 2.44e-17 2.45e-17 2.36e-17 2.45e-17 ... \cr
  ..$ : num [1:51, 1:1221] 1.63e-16 5.16e-17 5.16e-17 5.16e-17 5.32e-17 ... \cr
 $ Zevo      :List of 13 \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 1 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.000238 0.011094 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 1 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.000234 0.010203 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 1 ... \cr
  .. ..$ SMgas : num [1:51] 0.00 0.00 0.00 0.00 4.15e-05 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.00016 0.00892 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 1 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.000179 0.008337 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 1 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.000147 0.008561 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 0.996 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0.00424 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0.000121 0.00791 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 1 0.993 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0 0.00748 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0 0.00478 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 1 0.998 0.988 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0 0.00216 0.01172 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0 0 0 0 0.00206 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 1 0.999 0.994 0.982 ... \cr
  .. ..$ SMgas : num [1:51] 0 0 0.00143 0.00621 0.01825 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0.00 0.00 0.00 2.69e-09 1.10e-03 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 0.994 0.99 0.983 0.972 ... \cr
  .. ..$ SMgas : num [1:51] 0 0.00603 0.00999 0.0167 0.02793 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 1.06e-09 1.06e-09 1.06e-09 1.06e-09 5.21e-06 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 0.995 0.991 0.983 0.969 ... \cr
  .. ..$ SMgas : num [1:51] 0 0.00458 0.0091 0.01709 0.03085 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 7.40e-10 7.40e-10 7.40e-10 7.40e-10 2.46e-09 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 0.986 0.979 0.968 0.954 ... \cr
  .. ..$ SMgas : num [1:51] 0 0.0137 0.0209 0.0323 0.046 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0.00 0.00 0.00 0.00 3.45e-11 ... \cr
  ..$ :'data.frame':	51 obs. of  5 variables: \cr
  .. ..$ SMstar: num [1:51] 1 0.988 0.979 0.968 0.956 ... \cr
  .. ..$ SMgas : num [1:51] 0 0.0123 0.0209 0.0319 0.0444 ... \cr
  .. ..$ SMtot : num [1:51] 1 1 1 1 1 1 1 1 1 1 ... \cr
  .. ..$ SFR   : num [1:51] 1 0 0 0 0 0 0 0 0 0 ... \cr
  .. ..$ SMrem : num [1:51] 0.0 0.0 0.0 2.1e-10 2.1e-10 ... \cr
}
\details{
Everything has been computed for a Chabrier IMF out to 300 Msol binary variant of BPASS v2.2.1.
}
\source{
https://bpass.auckland.ac.nz/9.html
}
\references{
Stanway and Eldridge, 2018, MNRAS, 479, 75
}
\examples{
data(BPASS)

#Spectrum of our Sun today (Z=0.02, Age=5 Gyrs) [nearest to 4.6 Gyrs]:
plot(BPASS$Wave, BPASS$Zspec[[11]][38,], log='xy', type='l', xlim=c(2e3,4e4),
ylim=c(1e-8,1e-4), xlab=BPASS$Labels$Wavelab, ylab=BPASS$Labels$Lumlab)
}
\keyword{datasets}
