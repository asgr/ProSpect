\name{magABcalc}
\alias{magABcalc}
\alias{CGScalc}
\alias{Janskycalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Flux and Magnitude
}
\description{
Converts spectrum to a flux measurement in a filter (CGS or Jansky) or AB magnitude.
}
\usage{
magABcalc(wave, flux, filter = "r_VST")
CGScalc(wave, flux, filter='r_VST')
Janskycalc(wave, flux, filter='r_VST')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; observed wavelength in Angstroms. For convenience, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{flux}.
}
  \item{flux}{
Numeric vector; flux in erg/s / cm^2 / Angstrom specral units.
}
  \item{filter}{
Name of filter to use.
}
}
\details{
Depending on the format of the target photometry it may be more appropriate to work in AB mag, CGS or Jansky. As a reminder, Jansky is simply CGS*1e23, which scales things to nearer unity for typical astrophysical sources.
}
\value{
Numeric scalar; flux in target units (i.e. AB mag, CGS or Jansky).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{bandpass}}, \code{\link{getfilt}}, \code{\link{SFHp4}}
}
\examples{
fluxBC03=Lum2Flux(BC03lr$Wave, BC03lr$Zspec[[5]][161,]*1e10)

magABcalc(fluxBC03)
CGScalc(fluxBC03)
Janskycalc(fluxBC03)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ photometry }% use one of  RShowDoc("KEYWORDS")
