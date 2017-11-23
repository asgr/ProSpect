\name{photom}
\alias{photom}
\alias{Lum2Flux}
\alias{Lum2FluxFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral Processing for Photometry
}
\description{
These functions can manipulate restframe spectra in useful ways, e.g. placing it at a target redshift with the correct observed luminosity and band stretching.
}
\usage{
photom(wave, lum, filters = "all", z = 0.1, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = "planck")
Lum2Flux(wave, lum, z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR,
OmegaR = 0, w0 = -1, wprime = 0, ref='planck')
Lum2FluxFactor(z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR,
OmegaR = 0, w0 = -1, wprime = 0, ref='planck')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Numeric vector; the spectral wavelength in Angstroms.
}
  \item{lum}{
Numeric vector; the spectral luminosity in Lsun / Angstrom (BC03lr / BC03hr and EMILES format).
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters.
}
  \item{z}{
Numeric scalar; cosmological redshift, where z must be > -1.
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.3).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
}
  \item{OmegaR}{
Numeric scalar; Omega Radiation today (default is 0, but OmegaM/3400 is typical).
}
  \item{w0}{
The value of dark energy equation of state at z=0. See \code{\link{cosgrow}} for more details.
}
  \item{wprime}{
The evolution term that governs how the dark energy equation of state evolves with redshift. See \code{\link{cosgrow}} for more details.
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck=Planck15 and WMAP=WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
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
Aaron Robotham
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
function (wave, lum, filters = "all", z = 0.1, H0 = 100, OmegaM = 0.3, 
    OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, 
    ref = "planck") 
{
    if (!is.vector(wave)) {
        if (dim(wave)[2] == 2) {
            lum = wave[, 2]
            wave = wave[, 1]
        }
    }
    observedspec = Lum2Flux(wave = wave, lum = lum, z = z, H0 = H0, 
        OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, 
        wprime = wprime, ref = ref)
    if (filters[1] == "all") {
        filters = cenwave$filter
    }
    mag = {
    }
    for (i in filters) {
        mag = c(mag, magABcalc(observedspec, filter = i))
    }
    return(cbind(cenwave[match(filters, cenwave$filter), ], mag = mag))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
