\name{paranal_transmission}
\alias{paranal_transmission}
\docType{data}
\title{
EAZY Paranal Sky Transmission
}
\description{
Paranal sky transmission curve that comes with EAZY as downloaded on 20/11/2018.
}
\usage{data("paranal_transmission")}
\format{
  A data frame with 11514 observations on the following 2 variables.
  \describe{
    \item{\code{wave}}{a numeric vector, wavelength in Angstroms}
    \item{\code{response}}{a numeric vector, the relative response}
  }
}
\details{
Settings used:

Altitude :  65.38
Airmass  :   1.10

Season      : Entire Year
Time Period : Entire Night

PWV         : Seasonal Average
Monthly Averaged Solar Flux : 130.00 s.f.u
 
Components included in the radiance model:

Scattered Moon Light:
separation of Sun and Moon       :  90.00 degrees
separation of Moon and object    :  45.00 degrees
altitude of Moon above horizon   :  45.00 degrees
distance to Moon                 :   1.0000
 
Star Light
 
Zodiacal Light:
ecliptic longitude of object : 135.00 degrees
ecliptic latiude of object   :  90.00 degrees

Molecular Emission of Lower Atmosphere
Emission Lines of Upper Atmosphere
Airglow Continuum

Wavelength grid:
Range: 300.00 nm - 3000.00 nm
Fixed spectral resolution R=5000

LSF convolution kernel: none
}
\source{
The new web site to get such curves appears to be:

https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
}
\references{
Brammer, van Dokkum & Coppi, 2008, ApJ, 686, 1503
}
\examples{
data(paranal_transmission)
plot(paranal_transmission, type = 'l', xlab = 'Wave / Ang', ylab = 'Transmission')
}
\keyword{datasets}
