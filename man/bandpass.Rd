\name{bandpass}
\alias{bandpass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Band Pass
}
\description{
Put a spectrum through a filter response, returning either the modified spectrum or the integrated flux.
}
\usage{
bandpass(wave, flux, filter, lum = T)
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
Two column matrix; column one must be flux in Angstroms, column two must be a filter response. These are the type returned by 
\code{\link{getfilt} (see Examples).}
}
  \item{lum}{
Boolean scalar; should the integrated flux be returned (TRUE) or the spectral flux (FALSE).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Either returns the filtered spectrum (\option{lum}=FALSE) or the integrated flux (\option{lum}=TRUE).
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
function (wave, flux, filter, lum = T) 
{
    tempfunc = approxfun(x = filter[, 1], y = abs(filter[, 2]))
    tempremap = tempfunc(wave)
    tempremap[is.na(tempremap)] = 0
    if (lum) {
        return(sum(tempremap * wave * flux, na.rm = TRUE)/sum(tempremap * 
            wave, na.rm = TRUE))
    }
    else {
        return(tempremap * wave * flux/sum(tempremap * wave, 
            na.rm = TRUE))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
