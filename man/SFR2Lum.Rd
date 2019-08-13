\name{SFR2Lum}
\alias{SFR2Lum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SFR to H[alpha] Luminosity
}
\description{
A best effort conversion betwen SFR and H[alpha] line luminsity. Mostly for use with \code{\link{emissionLines}}.
}
\usage{
SFR2Lum(SFR = 1, lumscale = 21612724)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFR}{
Numeric scalar; star formation rate (Msol/yr).
}
  \item{lumscale}{
Numeric scalar; given by (1/7.9e-42) (Kennicutt, in erg/s) x 1e-7 (erg/s to Watt) x 2.61233e-27 (Watts to L[sol]) / 1.53 (Salpeter to Chabrier IMF conversion). Unless you have a very strong opinion on this scaling, the defaults is a reasonable best effort.
}
}
\details{
Best effort SFR to H[alpha] luminosity conversion function.
}
\value{
The H[alpha] line luminosity in Lsol.
}
\references{
Kennicutt et al, 1994, ApJ, 435, 22
Kennicutt et al, 1998, ARAA, 36, 189
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{emissionLines}}
}
\examples{
#Conversion of SFR=1 Msol/yr the H[alpha] line luminosity:

(1/7.9e-42)*1e-7*2.61233e-27/1.53
SFR2Lum()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ SFR }% use one of  RShowDoc("KEYWORDS")
