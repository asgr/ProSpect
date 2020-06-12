\name{Fritz}
\alias{Fritz}
\docType{data}
\title{
Fritz 2006 AGN Template Data
}
\description{
The full extended Fritz 2006 AGN model as per Feltre 2012. This is then interpolated by \code{\link{AGNinterp}} to create arbitrary AGN at non gridded points.
}
\usage{data("Fritz")}
\format{
  The format is:\cr
List of 9\cr
 $ Wave  : num [1:178] 10 11.2 12.6 14.1 15.8 ...\cr
 $ Aspec   : num [1:3, 1:4, 1:5, 1:8, 1:5, 1:10, 1:178] 0.0672 0.0698 0.0727 0.064 0.064 ...\cr
  ..- attr(*, "dimnames")=List of 7\cr
  .. ..$ : chr [1:3] "60" "100" "140"\cr
  .. ..$ : chr [1:4] "0" "2" "4" "6"\cr
  .. ..$ : chr [1:5] "-1" "-0.75" "-0.5" "-0.25" ...\cr
  .. ..$ : chr [1:8] "0.1" "0.3" "0.6" "1" ...\cr
  .. ..$ : chr [1:5] "10" "30" "60" "100" ...\cr
  .. ..$ : chr [1:10] "0" "10" "20" "30" ...\cr
  .. ..$ : NULL\cr
 $ params: chr [1:6] "ct" "al" "be" "ta" ...\cr
 $ ct    : num [1:3] 60 100 140\cr
 $ al    : num [1:4] 0 2 4 6\cr
 $ be    : num [1:5] -1 -0.75 -0.5 -0.25 0\cr
 $ ta    : num [1:8] 0.1 0.3 0.6 1 2 3 6 10\cr
 $ rm    : num [1:5] 10 30 60 100 150\cr
 $ an    : num [1:10] 0 10 20 30 40 50 60 70 80 90\cr
}
\details{
The below is taken directly from the included documentation:

A composition of power-laws, with different spectral indexes as a function of the wavelength, was used to model the primary source (i.e. accretion disk SED) SED taken from the file pwlbb05.dat.                                          

Standard, Galactic composition of dust, with graphite and astronomical silicates. The grain size distribution is a power law with index -3.5. Scattering and absorption coefficients taken from Laor and Draine (1993). Number of different species: 2 (Graphite and Silicate only).                      

The inner radius is set by the sublimation temperature of graphite grains: T = 1500 K. The dust sublimation radius is Rmin = 1.3*L46^0.5*T^(-2.8) pc [Barvainis (1987)] where L46=L/(1e46 erg/s), and T is the sublimation temperature. The dust density throughout the torus is expressed by the following: rho(r,theta)=A*r^(beta)exp(-gamma*cos(theta)), where A is a normalization costant, beta and gamma are, in the name of the model, be and al, respectively.

Each AGN model has been normalised so that the total luminosity fo the source AGN integrates to 1. Depending on the viewing angle the intergral of the observed AGN + torus will not in turn integrate to 1 (i.e. is is much less than this when viewed edge on).

\describe{
\item{Wave}{
  The wavelength in Angstroms [178 elements].
}
\item{Aspec}{
  7 dimensional array containing (in order) dimension of \option{ct} (opening angle of torus in degrees); \option{al} (gamma dust parameter); \option{be} (beta dust parameter); \option{ta} (optical depth tau); \option{rm} (outer to inner torus radius ratio); \option{an} (angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned), and the final column is the normalised luminosity per unit wavelength [178 elements].
}
\item{params}{
  Character vector; order list of parameters in \option{Aspec} model array.
}
\item{ct}{
  Numeric vector; available grid values of \option{ct}. Opening angle of torus in degrees.
}
\item{al}{
  Numeric vector; available grid values of \option{al}. Gamma dust parameter.
}
\item{be}{
  Numeric vector; available grid values of \option{be}. Beta dust parameter.
}
\item{ta}{
  Numeric vector; available grid values of \option{ta}. Optical depth tau.
}
\item{rm}{
  Numeric vector; available grid values of \option{rm}. Outer to inner torus radius ratio.
}
\item{an}{
  Numeric vector; available grid values of \option{an}. Angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned.
}
}
}
\source{
Downloaded the file ct20/ct40/ct60 from http://www.crya.unam.mx/gente/j.fritz/JFhp/AGN_models.html.
}
\references{
Feltre et al, 2012, MNRAS, 426, 120 \cr
Fritz et al, 2006, MNRAS, 366, 767
}
\examples{
data(Fritz)
str(Fritz)

plot(Fritz$Wave, Fritz$Aspec["100","4","-0.5","2","60","30",],
type='l', log='xy', xlab='Wave / Ang', ylab='Lum / L/A')
}
\keyword{datasets}
