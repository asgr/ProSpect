\name{magABcalc}
\alias{magABcalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Flux and Magnitude
}
\description{
Converts spectrum to measurement in filter (CGS or Jansky) or AB magnitude.
}
\usage{
magABcalc(wave, flux, filter = "r_VST", wavefac = 1e-10)
CGScalc(wave, flux, filter='r_VST', wavefac=1e-10)
Janskycalc(wave, flux, filter='r_VST', wavefac=1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Observed wavelength in Angstroms by default (see \option{wavefac}).
}
  \item{flux}{
Flux in usual erg/s / cm^2 / Angstrom specral units.
}
  \item{filter}{
Name of filter to use.
}
  \item{wavefac}{
Specified unit type converting 1m to chosen unit. Default of 1e-10 means Angstrom, whilst 1e-9 would be nm etc.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (wave, flux, filter = "r_VST", wavefac = 1e-10)
{
    c = 299792458
    if (!is.vector(wave)) {
        if (dim(wave)[2] == 2) {
            flux = wave[, 2]
            wave = wave[, 1]
        }
    }
    filter = getfilt(filter)
    fluxnu = (wavefac * flux * wave^2)/c
    totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter,
        lum = T)
    return(-2.5 * log10(totlumnu) - 48.6)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
