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
Either filter name explictly built in to \code{ProSpect}, or an unambiguous grep match to one of the available filters in \code{\link{EAZY_filters}}, e.g. 'cfht_mega_u' unabmigously matches to the full filter name of 'megaprime/cfht_mega_u_cfh9301.dat'. You can also conveniently specify your own tophat with the filter name syntax 'tophat_lo_hi' where 'lo' and 'hi' corresponds to the lower and upper limits of the filter tophat in Angstrom. \cr
For the full set now available see options at \code{\link{filters}}, where all listed filers are accessible to the \code{getfilt} function (with the removal of the preceding "filt_"). The current set available is: \cr
%%     ~~Describe \code{filter} here~~
\tabular{ll}{
  Name    \tab  Filter reference \cr

  % GALEX
  FUV_GALEX     \tab  GALEX FUV band \cr
  NUV_GALEX     \tab  GALEX NUV band \cr

  % SDSS
  u_SDSS  \tab  Sloan u band \cr
  g_SDSS  \tab  Sloan g band \cr
  r_SDSS  \tab  Sloan r band \cr
  i_SDSS  \tab  Sloan i band \cr
  z_SDSS  \tab  Sloan z band \cr

  % VST
  u_VST  \tab  VST u band \cr
  g_VST  \tab  VST g band \cr
  r_VST  \tab  VST r band \cr
  i_VST  \tab  VST i band \cr
  z_VST  \tab  VST z band \cr

  % LSST / Rubin
  u_LSST  \tab  Rubin LSST u band \cr
  g_LSST  \tab  Rubin LSST g band \cr
  r_LSST  \tab  Rubin LSST r band \cr
  i_LSST  \tab  Rubin LSST i band \cr
  z_LSST  \tab  Rubin LSST z band \cr
  y_LSST  \tab  Rubin LSST y band \cr

  % HSC
  g_HSC  \tab  HSC g band \cr
  r_HSC  \tab  HSC r band \cr
  i_HSC  \tab  HSC i band \cr
  z_HSC  \tab  HSC z band \cr
  Y_HSC  \tab  HSC Y band \cr

  % VISTA
  Z_VISTA  \tab  VISTA Z band \cr
  Y_VISTA  \tab  VISTA Y band \cr
  J_VISTA  \tab  VISTA J band \cr
  H_VISTA  \tab  VISTA H band \cr
  K_VISTA  \tab  VISTA K band \cr
  Ks_VISTA \tab  VISTA Ks band \cr

  % UKIRT
  Z_UKIRT  \tab  UKIRT Z band \cr
  Y_UKIRT  \tab  UKIRT Y band \cr
  J_UKIRT  \tab  UKIRT J band \cr
  H_UKIRT  \tab  UKIRT H band \cr
  K_UKIRT  \tab  UKIRT K band \cr

  % JWST (NIRCam + MIRI)
  F070W_JWST \tab JWST NIRCam F070W \cr
  F090W_JWST \tab JWST NIRCam F090W \cr
  F115W_JWST \tab JWST NIRCam F115W \cr
  F140M_JWST \tab JWST NIRCam F140M \cr
  F150W_JWST \tab JWST NIRCam F150W \cr
  F150W2_JWST \tab JWST NIRCam F150W2 \cr
  F162M_JWST \tab JWST NIRCam F162M \cr
  F164N_JWST \tab JWST NIRCam F164N \cr
  F182M_JWST \tab JWST NIRCam F182M \cr
  F187N_JWST \tab JWST NIRCam F187N \cr
  F200W_JWST \tab JWST NIRCam F200W \cr
  F210M_JWST \tab JWST NIRCam F210M \cr
  F212N_JWST \tab JWST NIRCam F212N \cr
  F250M_JWST \tab JWST NIRCam F250M \cr
  F277W_JWST \tab JWST NIRCam F277W \cr
  F300M_JWST \tab JWST NIRCam F300M \cr
  F322W2_JWST \tab JWST NIRCam F322W2 \cr
  F323N_JWST \tab JWST NIRCam F323N \cr
  F335M_JWST \tab JWST NIRCam F335M \cr
  F356W_JWST \tab JWST NIRCam F356W \cr
  F360M_JWST \tab JWST NIRCam F360M \cr
  F405N_JWST \tab JWST NIRCam F405N \cr
  F410M_JWST \tab JWST NIRCam F410M \cr
  F430M_JWST \tab JWST NIRCam F430M \cr
  F444W_JWST \tab JWST NIRCam F444W \cr
  F460M_JWST \tab JWST NIRCam F460M \cr
  F466N_JWST \tab JWST NIRCam F466N \cr
  F470N_JWST \tab JWST NIRCam F470N \cr
  F480M_JWST \tab JWST NIRCam F480M \cr
  F560W_JWST \tab JWST MIRI F560W \cr
  F770W_JWST \tab JWST MIRI F770W \cr
  F1000W_JWST \tab JWST MIRI F1000W \cr
  F1130W_JWST \tab JWST MIRI F1130W \cr
  F1280W_JWST \tab JWST MIRI F1280W \cr
  F1500W_JWST \tab JWST MIRI F1500W \cr
  F1800W_JWST \tab JWST MIRI F1800W \cr
  F2100W_JWST \tab JWST MIRI F2100W \cr
  F2550W_JWST \tab JWST MIRI F2550W \cr

  % Euclid
  VIS_Euclid   \tab Euclid VIS band \cr
  Y_Euclid     \tab Euclid Y band \cr
  J_Euclid     \tab Euclid J band \cr
  H_Euclid     \tab Euclid H band \cr
  Blue_Euclid  \tab Euclid Blue band \cr
  Red_Euclid   \tab Euclid Red band \cr

  % WFIRST / Roman
  R062_WFIRST  \tab Roman R062 band \cr
  Z087_WFIRST  \tab Roman Z087 band \cr
  Y106_WFIRST  \tab Roman Y106 band \cr
  J129_WFIRST  \tab Roman J129 band \cr
  W146_WFIRST  \tab Roman W146 band \cr
  H158_WFIRST  \tab Roman H158 band \cr
  F184_WFIRST  \tab Roman F184 band \cr

  % WISE
  W1_WISE \tab WISE W1 band \cr
  W2_WISE \tab WISE W2 band \cr
  W3_WISE \tab WISE W3 band \cr
  W4_WISE \tab WISE W4 band \cr

  % Spitzer
  I1_Spitzer \tab Spitzer IRAC I1 \cr
  I2_Spitzer \tab Spitzer IRAC I2 \cr
  I3_Spitzer \tab Spitzer IRAC I3 \cr
  I4_Spitzer \tab Spitzer IRAC I4 \cr
  M24_Spitzer \tab Spitzer MIPS 24 \cr
  M70_Spitzer \tab Spitzer MIPS 70 \cr
  M160_Spitzer \tab Spitzer MIPS 160 \cr

  % Herschel
  P70_Herschel \tab Herschel PACS 70 \cr
  P100_Herschel \tab Herschel PACS 100 \cr
  P160_Herschel \tab Herschel PACS 160 \cr
  S250_Herschel \tab Herschel SPIRE 250 \cr
  S350_Herschel \tab Herschel SPIRE 350 \cr
  S500_Herschel \tab Herschel SPIRE 500 \cr

  % JCMT
  S450_JCMT \tab JCMT SCUBA-2 450 \cr
  S850_JCMT \tab JCMT SCUBA-2 850 \cr

  % mm
  1mm_Aztec \tab AzTEC 1mm \cr
  2mm_Gizmo \tab GISMO 2mm \cr

  % ALMA
  Band9_ALMA \tab ALMA Band 9 \cr
  Band8_ALMA \tab ALMA Band 8 \cr
  Band7_ALMA \tab ALMA Band 7 \cr
  Band6_ALMA \tab ALMA Band 6 \cr
  Band5_ALMA \tab ALMA Band 5 \cr
  Band4_ALMA \tab ALMA Band 4 \cr
  Band3_ALMA \tab ALMA Band 3 \cr
  Band2_ALMA \tab ALMA Band 2 \cr
  Band1_ALMA \tab ALMA Band 1 \cr

  % VLA
  BandQ_VLA \tab VLA Q band \cr
  BandKa_VLA \tab VLA Ka band \cr
  BandK_VLA \tab VLA K band \cr
  BandKu_VLA \tab VLA Ku band \cr
  BandX_VLA \tab VLA X band \cr
  BandC_VLA \tab VLA C band \cr
  BandS_VLA \tab VLA S band \cr
  BandL_VLA \tab VLA L band \cr
  BandP_VLA \tab VLA P band \cr

  % GMRT
  Band5_GMRT \tab GMRT Band 5 \cr
  Band4_GMRT \tab GMRT Band 4 \cr
  Band3_GMRT \tab GMRT Band 3 \cr
  Band2_GMRT \tab GMRT Band 2 \cr
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
