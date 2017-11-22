\name{filters}
\alias{filt_g_SDSS}
\alias{filt_g_VST}
\alias{filt_H_VISTA}
\alias{filt_i_SDSS}
\alias{filt_i_VST}
\alias{filt_J_VISTA}
\alias{filt_K_VISTA}
\alias{filt_Ks_VISTA}
\alias{filt_NUV_GALEX}
\alias{filt_P100_Herschel}
\alias{filt_P160_Herschel}
\alias{filt_r_SDSS}
\alias{filt_r_VST}
\alias{filt_S250_Herschel}
\alias{filt_S350_Herschel}
\alias{filt_S450_Herschel}
\alias{filt_u_SDSS}
\alias{filt_u_VST}
\alias{filt_W1_WISE}
\alias{filt_W2_WISE}
\alias{filt_W3_WISE}
\alias{filt_W4_WISE}
\alias{filt_Y_VISTA}
\alias{filt_z_SDSS}
\alias{filt_Z_VISTA}
\alias{filt_z_VST}
\alias{filt}
\alias{filter}
\alias{filters}
\docType{data}
\title{
Band Filters
}
\description{
Various common filter response functions, appropriate for using in calculating magitudes etc.
}
\usage{
data("filt_FUV_GALEX")
data("filt_g_SDSS")
data("filt_g_VST")
data("filt_H_VISTA")
data("filt_i_SDSS")
data("filt_i_VST")
data("filt_J_VISTA")
data("filt_K_VISTA")
data("filt_Ks_VISTA")
data("filt_NUV_GALEX")
data("filt_P100_Herschel")
data("filt_P160_Herschel")
data("filt_r_SDSS")
data("filt_r_VST")
data("filt_S250_Herschel")
data("filt_S350_Herschel")
data("filt_S450_Herschel")
data("filt_u_SDSS")
data("filt_u_VST")
data("filt_W1_WISE")
data("filt_W2_WISE")
data("filt_W3_WISE")
data("filt_W4_WISE")
data("filt_Y_VISTA")
data("filt_z_SDSS")
data("filt_Z_VISTA")
data("filt_z_VST")
}
\format{
  A number of data frames with the following 2 variables.
  \describe{
    \item{\code{wave}}{A numeric vector}
    \item{\code{response}}{A numeric vector}
  }
}
\details{
Note the filter for VISTA is actually Ks (for K-short) which makes it approximately similar to 2MASS. This differs to WFCAM (e.g. UKIDSS) which used a full K filter. Both K_VISTA and Ks_VISTA are provided to prevent silly name convention mistakes, but these response curves are in fact identical (i.e. Ks).
}
\source{
Mostly the relevant websites for various facilities.
}
\examples{
colvec=rev(rainbow(9, end=2/3))
plot(filt_u_VST, type='l', col=colvec[1], xlim=c(3e3,2.5e4), ylim=c(0,1))
lines(filt_g_VST, type='l', col=colvec[2])
lines(filt_r_VST, type='l', col=colvec[3])
lines(filt_i_VST, type='l', col=colvec[4])
lines(filt_Z_VISTA, type='l', col=colvec[5])
lines(filt_Y_VISTA, type='l', col=colvec[6])
lines(filt_J_VISTA, type='l', col=colvec[7])
lines(filt_H_VISTA, type='l', col=colvec[8])
lines(filt_K_VISTA, type='l', col=colvec[9])
}
\keyword{datasets}
