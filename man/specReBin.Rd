\name{specReBin}
\alias{specReBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral Re-Binner
}
\description{
It is often useful to rebin spectra to either log binning or linear binning. This function does just that.
}
\usage{
specReBin(wave, flux, wavegrid = NULL, bin = NULL, binfunc = median, interp = 'approx',
  logbin = TRUE, rough = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; the input spectral wavelength in any linear units. For convenience, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{lum}.
}
  \item{flux}{
Numeric vector; the input flux density in any linear units.
}
  \item{wavegrid}{
Numeric vector; the target spectral wavelength grid in any linear units (must be same as \option{wave}).
}
  \item{bin}{
Numeric scalar; the bin size to use in either liner units (same as \option{wave}) or log units (i.e. dex). If not provided this will be automatically computed.
}
  \item{binfunc}{
Function; if \option{bin} is not provived this function will be run on either the linear or logged wavelength (as per \option{logbin}) to choose a reasonable constant bin size.
}
  \item{interp}{
Character scalar; The interpolator to use. Either 'approx' (see \code{\link{approx}}) or 'spline' (see \code{\link{spline}}). The former is better behaved with discrete spectrum features, but the latter might work better when the spectrum is sufficiently smooth.
}
  \item{logbin}{
Logical; should the automatic wavelength grid be evenly log spaced? If TRUE it will be evenly log spaced, if FALSE it will be evenly linearly spaced.
}
  \item{rough}{
Logical; should rough interpolation be used to rebin the spectrum (default is TRUE) or should full be pixel integration be done (FALSE). The latter is usually only necessary when degrading the spectrum to coarser resolution because it ensures strong emission line features are correctly integrated over during re-binning. Note the full integration is very slow (takes about a second for a typical spectrum), so you only want to use this option when critically necessary.
}
  \item{\dots}{
Extra arguments to pass to either \code{\link{approx}} or \code{\link{spline}} as per \option{interp}.
}
}
\details{
It is often useful to rebin spectra to either log binning or linear binning. This function does just that.

Note you need to be careful the resolution is high enough to not lose flux. This function just interpolates the flux density, it does not integrate over bins and reditribute flux. To do a more accurate job you would rebin at very high resolution and then integrate that spectrum down into coarser bins.
}
\value{
Data.frame; first column is wavelength and second column is observed flux density.
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{approx}}, \code{\link{spline}}, \code{\link{addspec}}
}
\examples{
library(magicaxis)
tempflux = ProSpectSED(sparse=1, emission=TRUE)$FinalFlux
tempflux = tempflux[tempflux$wave>2e3 & tempflux$wave<2e4,]

magplot(tempflux, type='l', xlim=c(7e3,8e3))
lines(specReBin(tempflux), col='red') #works pretty well automatically
lines(specReBin(tempflux, bin=1e-3), col='green') #mmisses lines
lines(specReBin(tempflux, bin=1e-5), col='blue') #looks very close to the original

# we can also compare rough and accuracte integration schemes:

magplot(tempflux, type='l', xlim=c(7e3,8e3))
lines(specReBin(tempflux), col='red') #rough is on by default
lines(specReBin(tempflux, rough=FALSE), col='blue') #accurate not any better (and slow)

}
\concept{ rebin }
