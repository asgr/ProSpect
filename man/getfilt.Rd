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
Either filter name explictly built in to \code{ProSpect}, or an unambiguous grep match to one of the available filters in \code{\link{EAZY_filters}}, e.g. 'cfht_mega_u' unabmigously matches to the full filter name of 'megaprime/cfht_mega_u_cfh9301.dat'. \cr
For the full set now available see options at \code{\link{filters}}, where all data with are accessible to the \code{getfilt} function (with the removal of the preceding "filt_"). The original set available was: \cr
%%     ~~Describe \code{filter} here~~
\tabular{ll}{
  Name    \tab  Filter reference \cr
  FUV     \tab  GALEX FUV band \cr
  NUV     \tab  GALEX NUV band \cr
  u  \tab  Sloan u band \cr
  g  \tab  Sloan g band \cr
  r  \tab  Sloan r band \cr
  i  \tab  Sloan i band \cr
  z  \tab  Sloan z band \cr
  u_VST  \tab  VST u band \cr
  g_VST  \tab  VST g band \cr
  r_VST  \tab  VST r band \cr
  i_VST  \tab  VST i band \cr
  z_VST  \tab  VST z band \cr
  g_HSC  \tab  HSC g band \cr
  r_HSC  \tab  HSC r band \cr
  i_HSC  \tab  HSC i band \cr
  z_HSC  \tab  HSC z band \cr
  Y_HSC  \tab  HSC Y band \cr
  Z  \tab  VISTA Z band \cr
  Y  \tab  VISTA Y band \cr
  J  \tab  VISTA J band \cr
  H  \tab  VISTA H band \cr
  K or Ks  \tab  VISTA Ks band \cr
  Z_UKIRT  \tab  UKIRT Z band \cr
  Y_UKIRT  \tab  UKIRT Y band \cr
  J_UKIRT  \tab  UKIRT J band \cr
  H_UKIRT  \tab  UKIRT H band \cr
  K_UKIRT  \tab  UKIRT K band \cr
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
  Band9   \tab  ALMA Band 9 \cr
  Band8   \tab  ALMA Band 8 \cr
  Band7   \tab  ALMA Band 7 \cr
  Band6   \tab  ALMA Band 6 \cr
  Band5   \tab  ALMA Band 5 \cr
  Band4   \tab  ALMA Band 4 \cr
}
}
}
\details{
Returns a target named filter from the internal database of filters. Filters returned are of the type to be used when calculating AB magnitudes in \code{\link{bandpass}}.
\cr
Note LSST filters can also be accessed with "Rubin" (the preferred new shorthand name for the observatory).
}
\value{
Two column matrix. First column is the wavelength in Angstrom and second column is the filter response.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{bandpass}}, \code{\link{magABcalc}}, \code{\link{SFHfunc}}
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
\concept{ filter }% use one of  RShowDoc("KEYWORDS")
