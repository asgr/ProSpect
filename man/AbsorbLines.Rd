\name{AbsorbLines}
\alias{AbsorbLines}
\docType{data}
\title{
Galaxy Absorption Lines
}
\description{
Significant Optical and NIR Galaxy Absorption Lines. This covers the famous A-K Fraunhofer lines, as well as other hydrogen and oxygen elemental lines.
}
\usage{data("AbsorbLines")}
\format{
  A data frame with 47 observations on the following 5 variables.
  \describe{
    \item{\code{wave}}{Numeric; wavelength (Ang)}
    \item{\code{width}}{Numeric; feature width (Ang)}
    \item{\code{element}}{Character; name of atomic element.}
    \item{\code{state}}{Character; atomic state (neutral line = I, single ionised = II, double ionised = III etc).}
    \item{\code{fraunhofer}}{Character; Fraunhofer label if relevant. These are only labelled C-K (A and B are not within the range listed here, and are terrestrial oxygen anyway).}
  }
}
\details{
The list of absorption features listed here is not necessarily comprehensive, but it covers many of the common lines seen in galaxy spectra. Note the G-band (sometimes stylised as the g-band) is not a single feature, rather a complex of methane (CH), iron (Fe) and hydrogen (H-delta). Only very high resolution data tends to separate out the various components. The dominated calcium H-band is also shared with the hydrogen feature (H-episilon).
}
\source{
Various online material (pretty standard list).
}
\examples{
data(AbsorbLines)

# The famous Calcium HK doublet:
AbsorbLines[AbsorbLines$fraunhofer == 'H',]
AbsorbLines[AbsorbLines$fraunhofer == 'K',]
}
\keyword{datasets}
