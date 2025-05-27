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

Altitude :  65.38 \cr
Airmass  :   1.10

Season      : Entire Year \cr
Time Period : Entire Night

PWV         : Seasonal Average \cr
Monthly Averaged Solar Flux : 130.00 s.f.u

Components included in the radiance model:

Scattered Moon Light:
separation of Sun and Moon       :  90.00 degrees \cr
separation of Moon and object    :  45.00 degrees \cr
altitude of Moon above horizon   :  45.00 degrees \cr
distance to Moon                 :   1.0000

Star Light

Zodiacal Light:
ecliptic longitude of object : 135.00 degrees \cr
ecliptic latiude of object   :  90.00 degrees

Molecular Emission of Lower Atmosphere \cr
Emission Lines of Upper Atmosphere \cr
Airglow Continuum

Wavelength grid: \cr
Range: 300.00 nm - 3000.00 nm \cr
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

# Use on a filter:
paranal_func=approxfun(paranal_transmission)

SDSS_u=EAZY_filters$filters[[73]]
plot(SDSS_u, type='l')

SDSS_u_sky=cbind(wave=SDSS_u[,1],response=SDSS_u[,2]*paranal_func(SDSS_u[,1]))
lines(SDSS_u_sky, lty=2)
}
\keyword{datasets}
