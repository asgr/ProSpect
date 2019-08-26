\name{addspec}
\alias{addspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add Mismatched Spectra Together
}
\description{
Simple function that interpolates two non wavelength alligned spectra onto a common interpolated scheme,
}
\usage{
addspec(wave1, flux1, wave2, flux2, extrap = "constant", waveout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave1}{
Numeric vector; wavelength 1 in in any units (but must be the same as \option{wave2}).
}
  \item{flux1}{
Numeric vector; flux 1 in any units (but must be the same as \option{flux2}).
}
  \item{wave2}{
Numeric vector; wavelength 2 in in any units (but must be the same as \option{wave1}).
}
  \item{flux2}{
Numeric vector; flux 2 in any units (but must be the same as \option{flux1}).
}
  \item{extrap}{
If 'constant' then values outside the range inherit the nearest. If a numeric scalar (e.g. 0) then values outside the range inherit this value. 
}
  \item{waveout}{
Numeric vector; desired output wavelength grid to use (but must be the same as \option{wave1}/\option{wave2}). If not provided then a sorted vector of \option{wave1} and \option{wave2} is used instead.
}
}
\details{
Internally it uses the R \code{\link{approxfun}}, so if more sophisticated interpolation is necessary then other routes might be required.
}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{approxfun}}
}
\examples{
plot(BC03lr$Wave, BC03lr$Zspec[[5]][100,], log='xy', type='l', xlim=c(1e3,1e7),
ylim=c(1e-4,1), xlab=BC03lr$Labels$Wavelab, ylab=BC03lr$Labels$Lumlab, col='grey')
lines(Dale_Msol$Wave, Dale_Msol$Aspec[[1]][10,], col='grey')
add=addspec(BC03lr$Wave, BC03lr$Zspec[[5]][100,],Dale_Msol$Wave,
Dale_Msol$Aspec[[1]][10,])
lines(add)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
