\name{getfilt}
\alias{getfilt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Requested Filter Curve
}
\description{
Returns a target named filter from the internal database of filters.
}
\usage{
getfilt(filter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
%%     ~~Describe \code{filter} here~~
\tabular{ll}{
  Name    \tab  Filter reference \cr
  FUV     \tab  GALEX FUV band \cr
  NUV     \tab  GALEX NUV band \cr
  u_SDSS  \tab  Sloan u band \cr
  g_SDSS  \tab  Sloan g band \cr
  r_SDSS  \tab  Sloan r band \cr
  i_SDSS  \tab  Sloan i band \cr
  z_SDSS  \tab  Sloan z band \cr
  u_VST  \tab  VST u band \cr
  g_VST  \tab  VST g band \cr
  r_VST  \tab  VST r band \cr
  i_VST  \tab  VST i band \cr
  z_VST  \tab  VST z band \cr
  Z_VISTA  \tab  VISTA Z band \cr
  Y_VISTA  \tab  VISTA Y band \cr
  J_VISTA  \tab  VISTA J band \cr
  H_VISTA  \tab  VISTA H band \cr
  K_VISTA or Ks_VISTA  \tab  VISTA Ks band \cr
  W1      \tab  Wise W1 band  \cr
  W2      \tab  Wise W2 band  \cr
  W3      \tab  Wise W3 band  \cr
  W4      \tab  Wise W4 band  \cr
  I1      \tab  Spitzer I1 band  \cr
  I2      \tab  Spitzer I2 band  \cr
  I3      \tab  Spitzer I3 band  \cr
  I4      \tab  Spitzer I4 band  \cr
  M24     \tab  Spitzer M24 band  \cr
  M70     \tab  Spitzer M70 band  \cr
  M160    \tab  Spitzer M160 band  \cr
  P70     \tab  Herschel PACS 70 band  \cr
  P100    \tab  Herschel PACS 100 band  \cr
  P160    \tab  Herschel PACS 160 band  \cr
  S250    \tab  Herschel Spire 250 band  \cr
  S350    \tab  Herschel Spire 350 band  \cr
  S500    \tab  Herschel Spire 500 band  \cr
  S450    \tab  JCMT SCUBA-2 450 band  \cr
  S850    \tab  JCMT SCUBA-2 850 band  \cr
  1mm     \tab  Aztec 1mm band \cr
  2mm     \tab  Gizmo 2mm band \cr
}
}
}
\details{
Returns a target named filter from the internal database of filters. Filters returned are of the type to be used when calculating AB magnitudes in \code{\link{bandpass}}.
}
\value{
Two column matrix. First column is the wavelength in Angstrom and second column is the filter response.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{bandpass}}, \code{\link{magABcalc}}, \code{\link{SFHp4}}
}
\examples{
colvec=rev(rainbow(9, end=2/3))
plot(getfilt('u_VST'), type='l', col=colvec[1], xlim=c(3e3,2.5e4), ylim=c(0,1),
xlab = 'Wave / Ang', ylab = 'Response')
lines(getfilt('g_VST'), type='l', col=colvec[2])
lines(getfilt('r_VST'), type='l', col=colvec[3])
lines(getfilt('i_VST'), type='l', col=colvec[4])
lines(getfilt('Z_VISTA'), type='l', col=colvec[5])
lines(getfilt('Y_VISTA'), type='l', col=colvec[6])
lines(getfilt('J_VISTA'), type='l', col=colvec[7])
lines(getfilt('H_VISTA'), type='l', col=colvec[8])
lines(getfilt('K_VISTA'), type='l', col=colvec[9])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
