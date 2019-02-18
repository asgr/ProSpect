\name{ProSpect-package}
\alias{ProSpect-package}
\alias{ProSpect}
\docType{package}
\title{
\packageTitle{ProSpect}
}
\description{
\packageDescription{ProSpect}
}
\details{
\tabular{ll}{
Package: \tab ProSpect\cr
Type: \tab Package\cr
Version: \tab 0.6.0\cr
Date: \tab 2019-02-18\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 2.13)\cr
Imports:  \tab celestial, NISTunits, pracma\cr
}
}

\author{
\packageAuthor{ProSpect}

Maintainer: \packageMaintainer{ProSpect}
}
\references{
Bertelli et al., 1994, A\&AS, 106, 275
Bruzual and Charlot, 2003, MNRAS, 344, 1000
Chabrier, 2003, PASP, 115, 763
Dale et al, 2014, ApJ, 784, 11
Vazekis et al, 2016, MNRAS, 463, 3409
}
\keyword{SED}

\examples{
SFHdemo_dust=SFHp4()
plot(SFHdemo_dust$out[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)
}
