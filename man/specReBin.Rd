\name{specReBin}
\alias{specReBin}
\alias{speclibReBin}
\alias{speclibReGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SSP Spectral Re-Binner and SSP Re-Gridder Functions
}
\description{
It is often useful to rebin spectra to either log binning or linear binning. This function does just that. Also allows an entire SSP library to be rebinned (useful is using a very high resolution SSP for broad-band work etc).
}
\usage{
specReBin(wave, flux, wavegrid = NULL, invar = NULL, bin = NULL, binfunc = median,
  interp = 'approx',  logbin = TRUE, rough = FALSE, ...)

speclibReBin(speclib, wavegrid = NULL, bin = NULL, binfunc = median, interp = 'approx',
  logbin = TRUE, rough = FALSE, cores = 1L, check = FALSE, ...)

speclibReGrid(speclib, logAge_steps = NULL, logZ_steps = NULL, Zsol = 0.02, cores = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; the input spectral wavelength in any linear units. For convenience, if a two column matrix or data frame is supplied then the first column is taken to be \option{wave} and the second is taken to be \option{lum}.
}
  \item{flux}{
Numeric vector; the input flux density in any linear units.
}
  \item{speclib}{
List; the spectral library (e.g. BC03hr).
}
  \item{wavegrid}{
Numeric vector; the target spectral wavelength grid in any linear units (must be same as \option{wave}). If this is missing then it will be computed automatically using \option{bin}. See also \option{binfunc}.
}
  \item{invar}{
Numeric vector; inverse variance weighting corresponding to the input \option{flux} (i.e. if provided this needs to be the same length). If left as NULL then this will be ignored, and weighting free re-binning will be used instead.
}
  \item{bin}{
Numeric scalar; the bin size to use in either liner units (same as \option{wave}) or log units (i.e. dex). If not provided this will be automatically computed. If this is provided and \option{wavegrid} is missing, then internally \option{wavegrid} is computed to be \code{seq(floor(wave), ceiling(wave), by = bin)}. See also \option{binfunc}.
}
  \item{binfunc}{
Function; if \option{bin} is not provived this function will be run on either the linear or logged wavelength (as per \option{logbin}) to choose a reasonable constant bin size.
}
  \item{interp}{
Character scalar; The interpolator to use. Either 'approx' (see \code{\link{approx}}) or 'spline' (see \code{\link{spline}}). The former is better behaved with discrete spectrum features, but the latter might work better when the spectrum is sufficiently smooth.
}
  \item{logbin}{
Logical; should the automatic wavelength grid be evenly log spaced? If TRUE it will be evenly log spaced, if FALSE it will be evenly linearly spaced. Say you have a wavelength 10 and 20, a linear bin would be at 15 and a log bin would be at 14.1.
}
  \item{rough}{
Logical; should rough interpolation be used to rebin the spectrum (TRUE) or should full pixel integration be done (FALSE, default). As of \code{ProSpect} v1.4.1 \option{rough} = FALSE is just as fast (if not faster) than the coarse interpolation method. For spectral fitting you certainly might want the more accurate integration, but it can cause discontinuities when moving from very coarse model data to quite fine resolution observed data.
}
  \item{cores}{
Integer scalar; the number of cores to use.
}
  \item{check}{
Logical; should checks be run on the spectral library being returned? Uses \code{speclib_check}. If checks fail the SSP is still returned, and can then be checked and potentially fixed by the user.
}
  \item{logAge_steps}{
Numeric vector; if provided it is the log10(Age/yr) values that you want the output SSP to be interpolated at. If this argument is NULL we use the sorted unique log10(Age) values from the \option{speclib} provided (the default). This means you can generate SSPs with fairly arbitrary age grids, with the caveat the \option{logAge_steps} must be fully inside the logAge extremes in the provided \option{speclib} (we can only interpolate, not extrapolate). If length 1 and \option{logZ_steps} is also length 1 then just the interpolated spectrum is returned (see Value below).
}
  \item{logZ_steps}{
Numeric vector; if provided it is the log10(Z/Zsol) values that you want the output SSP to be interpolated at. If this argument is NULL we use the sorted unique logZ values from the \option{speclib} provided (the default). This means you can generate SSPs with fairly arbitrary metallicity grids, with the caveat the \option{logZ_steps} must be fully inside the logZ extremes in the provided \option{speclib} (we can only interpolate, not extrapolate). If length 1 and \option{logAge_steps} is also length 1 then just the interpolated spectrum is returned (see Value below).
}
  \item{Zsol}{
Numeric scalar; value of the Sun's solar metallicity (to convert logZ to Z).
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
\code{specReBin} Data.frame; first column is wavelength [wave] and second column is observed flux density [flux].

\code{speclibReBin} List; the re-binnned SSP spectral library.

\code{speclibReGrid} List; the re-gridded (to different Age and Z values) SSP spectral library. If the input \option{logAge_steps} and \option{logZ_steps} are both length 1 then as a convenience feature just the interpolated spectrum is returned in a data.frame with columns 'wave' (Ang) and 'lum' (Lsol/Ang).
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{approx}}, \code{\link{spline}}, \code{\link{addspec}}
}
\examples{
library(magicaxis)
tempflux = ProSpectSED(sparse=1, emission=TRUE, stellpop = "BC03hr")$FinalFlux
tempflux = tempflux[tempflux$wave>2e3 & tempflux$wave<2e4,]

magplot(tempflux, type='l', xlim=c(7e3,8e3))
lines(specReBin(tempflux, rough=TRUE), col='red') #works pretty well automatically
lines(specReBin(tempflux, bin=1e-3, rough=TRUE), col='green') #misses lines
lines(specReBin(tempflux, bin=1e-5, rough=TRUE), col='blue') #very close to the original

# We can also compare rough and accuracte integration schemes:

magplot(tempflux, type='l', xlim=c(7e3,8e3))
lines(specReBin(tempflux, bin=1e-3, rough=TRUE), col='red') #rough and fast
lines(specReBin(tempflux, bin=1e-3, rough=FALSE), col='blue') #conserves flux

# And now a whole spectral library:

data(BC03hr)
data(BC03lr)
speclib_rebin = speclibReBin(BC03hr, BC03lr$Wave, rough=FALSE)

# Note the ProSpect library re-bin seems to be more accurate than the official BC03 one!

magplot(BC03hr$Wave, BC03hr$Zspec[[1]][100,], type='l', log='xy', xlim=c(1e3,1e4))
lines(BC03lr$Wave, BC03lr$Zspec[[1]][100,], col='red')
lines(speclib_rebin$Wave, speclib_rebin$Zspec[[1]][100,], col='blue')
}
\concept{ rebin }
