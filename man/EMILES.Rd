\name{EMILES}
\alias{EMILES}
\docType{data}
\title{
E-MILES High Resolution Models
}
\description{
This is a reformatted collection of all the key informtion for the high resolution E-MILES models. This gives enough information to generate generic spectra for a target metallicity and star formation history.
}
\usage{data("EMILES")}

\format{
The highest level list structure contains the following elements:

\describe{
  \item{Z}{Numeric vector of length 12 describing the availability and location of different metallicity stellar populations.}
  \item{Age}{Numeric vector of length 53 giving the mid point of the age in years for the stellar populations.}
  \item{AgeBins}{Numeric vector of length 54 giving the age bin extremes in years of the stellar populations.}
  \item{AgeWeights}{Numeric vector of length 53 giving the total time in years for the stellar population bins.}
  \item{Wave}{Numeric vector of length 53689 giving the spectral wavelength in Angstroms.}
  \item{Labels}{Character vector of length 5 giving handy axis labels.}
  \item{Zpsec}{List containing 12 (for each metallicity) 53x53689 matrices that give the full spectral information.}
  \item{Zevo}{List containing 12 (for each metallicity) 53x5 data frames that contain the evolution information for the stellar material (stellar mass in stars [1] / gas [2] / total [3]; the star formation rate [4] and the stellar mass in remnants [5])}
}

The full format is:\cr
\cr
List of 8\cr
 $ Z         : num [1:12] 0.0001 0.0003 0.0006 0.001 0.002 0.004 0.008 0.01 0.0198 0.024 ...\cr
 $ Age       : num [1:53] 3.0e+07 4.0e+07 5.0e+07 6.0e+07 7.0e+07 8.0e+07 9.0e+07 1.0e+08 1.5e+08 2.0e+08 ...\cr
 $ AgeBins   : num [1:54] 0 0.005 0.015 0.025 0.035 0.045 0.055 0.065 0.095 0.145 ...\cr
 $ AgeWeights: num [1:53] 0.005 0.01 0.01 0.01 0.01 0.01 0.01 0.03 0.05 0.05 ...\cr
 $ Wave      : num [1:53689] 1680 1681 1682 1683 1684 ...\cr
 $ Labels    :List of 5\cr
  ..$ Zlab      : chr "Metallicity"\cr
  ..$ Agelab    : chr "Time since ZAM / Yrs"\cr
  ..$ Wavelab   : chr "Wavelength / Ang"\cr
  ..$ Lumlab   : chr "Lsol / Ang (for 1 Msol SF)"\cr
  ..$ LumAgelab: chr "Lsol / Ang (for 1 Msol/Yr SFR)"\cr
 $ Zspec     :List of 12\cr
  ..$ : num [1:53, 1:53689] 0.0265 0.0193 0.0151 0.0127 0.0111 ...\cr
  ..$ : num [1:53, 1:53689] 0.0259 0.0184 0.0144 0.0122 0.0107 ...\cr
  ..$ : num [1:53, 1:53689] 0.0249 0.0174 0.0139 0.0119 0.0104 ...\cr
  ..$ : num [1:53, 1:53689] 0.02293 0.01687 0.01337 0.01109 0.00953 ...\cr
  ..$ : num [1:53, 1:53689] 0.01999 0.01508 0.01214 0.00979 0.00817 ...\cr
  ..$ : num [1:53, 1:53689] 0.01705 0.01247 0.00976 0.00793 0.0067 ...\cr
  ..$ : num [1:53, 1:53689] 0.01576 0.01115 0.00847 0.00667 0.00535 ...\cr
  ..$ : num [1:53, 1:53689] 0.01601 0.01131 0.00837 0.00643 0.00521 ...\cr
  ..$ : num [1:53, 1:53689] 0.01378 0.00923 0.00674 0.00511 0.00402 ...\cr
  ..$ : num [1:53, 1:53689] 0.01276 0.0084 0.00619 0.00459 0.00377 ...\cr
  ..$ : num [1:53, 1:53689] 0.01187 0.00788 0.0057 0.00416 0.00344 ...\cr
  ..$ : num [1:53, 1:53689] 0.01073 0.00718 0.00481 0.00381 0.00329 ...\cr
 $ Zevo      :List of 12\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.786 0.767 0.752 0.74 0.731 ...\cr
  .. ..$ SMgas : num [1:53] 0.169 0.185 0.197 0.206 0.214 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0447 0.0479 0.0508 0.0533 0.0556 0.0576 0.0595 0.0614 0.0688 0.0736 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.766 0.751 0.739 0.729 ...\cr
  .. ..$ SMgas : num [1:53] 0.17 0.186 0.198 0.207 0.215 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0448 0.0481 0.051 0.0536 0.0559 0.058 0.06 0.0618 0.0685 0.0731 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.765 0.75 0.738 0.728 ...\cr
  .. ..$ SMgas : num [1:53] 0.171 0.187 0.199 0.208 0.216 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0449 0.0482 0.0512 0.0538 0.0561 0.0582 0.0603 0.062 0.0682 0.0729 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.784 0.765 0.75 0.738 0.728 ...\cr
  .. ..$ SMgas : num [1:53] 0.171 0.187 0.199 0.208 0.216 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0449 0.0483 0.0512 0.0538 0.0561 0.0583 0.0603 0.062 0.068 0.0725 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.784 0.765 0.75 0.738 0.728 ...\cr
  .. ..$ SMgas : num [1:53] 0.171 0.187 0.199 0.208 0.215 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0449 0.0482 0.0511 0.0537 0.0561 0.0582 0.0602 0.0619 0.0673 0.0718 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.765 0.751 0.739 0.729 ...\cr
  .. ..$ SMgas : num [1:53] 0.17 0.186 0.198 0.207 0.215 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0449 0.0482 0.051 0.0535 0.0558 0.0579 0.0599 0.0616 0.067 0.0711 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.766 0.752 0.741 0.731 ...\cr
  .. ..$ SMgas : num [1:53] 0.17 0.186 0.197 0.206 0.213 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0448 0.048 0.0508 0.0532 0.0554 0.0574 0.0593 0.0611 0.0665 0.0708 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.767 0.753 0.742 0.732 ...\cr
  .. ..$ SMgas : num [1:53] 0.17 0.185 0.197 0.205 0.213 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0448 0.0479 0.0506 0.053 0.0552 0.0572 0.059 0.0608 0.066 0.0701 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.785 0.767 0.754 0.743 0.734 ...\cr
  .. ..$ SMgas : num [1:53] 0.17 0.185 0.196 0.204 0.211 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0448 0.0478 0.0505 0.0527 0.0547 0.0566 0.0583 0.0599 0.065 0.0688 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.784 0.767 0.753 0.743 0.734 ...\cr
  .. ..$ SMgas : num [1:53] 0.171 0.185 0.196 0.204 0.211 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0449 0.0479 0.0505 0.0527 0.0547 0.0565 0.0582 0.0598 0.0649 0.0686 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.783 0.766 0.753 0.743 0.734 ...\cr
  .. ..$ SMgas : num [1:53] 0.172 0.186 0.196 0.204 0.211 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0451 0.048 0.0506 0.0527 0.0547 0.0565 0.0581 0.0597 0.0647 0.0683 ...\cr
  ..$ :'data.frame':	53 obs. of  5 variables:\cr
  .. ..$ SMstar: num [1:53] 0.781 0.764 0.752 0.742 0.733 ...\cr
  .. ..$ SMgas : num [1:53] 0.173 0.187 0.198 0.205 0.212 ...\cr
  .. ..$ SMtot : num [1:53] 1 1 1 1 1 1 1 1 1 1 ...\cr
  .. ..$ SFR   : num [1:53] 1 0 0 0 0 0 0 0 0 0 ...\cr
  .. ..$ SMrem : num [1:53] 0.0455 0.0484 0.0509 0.053 0.055 0.0567 0.0584 0.06 0.0649 0.0685 ...\cr
}
\details{
There are 12 metallicities (as given in \option{Z}), 53 age steps per population and 53689 spectral fluxes. Everything has been computed for a Chabrier IMF and the BaSTI isochrones.
}
\source{
http://miles.iac.es/
}
\references{
EMILES: Vazdekis et al, 2016, MNRAS, 463, 3409 \cr
Padova 1994 isochrones: Bertelli et al., 1994, A&AS, 106, 275 \cr
Chabrier IMF: Chabrier, 2003, PASP, 115, 763
}
\examples{
#The spectrum of our Sun today (Z=0.02, Age=4.5 Gyrs [nearest to 4.6 Gyrs]):
plot(EMILES$Wave, EMILES$Zspec[[9]][34,], log='xy', type='l', xlim=c(2e3,4e4),
ylim=c(1e-8,1e-4), xlab=EMILES$Labels$Wavelab, ylab=EMILES$Labels$Lumlab)
}
\keyword{datasets}
