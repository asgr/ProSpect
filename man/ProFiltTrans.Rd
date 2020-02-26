\name{ProFiltTrans}
\alias{ProFiltTrans}
\alias{ProFiltTrans_Ran}
\alias{ProFiltTrans_Shark}
\docType{data}
\title{
Pre-computed Magnitudes for Random SFHs
}
\description{
A list of pre-computed magnitudes over a range of lookback time (separated in Gyr intervals) for randomly sampled SFHs using \code{\link{massfunc_b5}} and Uniform sampling (from 0 to 1 Msol/yr) for each mass bin. These outputs can be used to compute filter conversions without the added cost and complexity of creating \code{\link{ProSpectSED}} SEDs.
}
\usage{
data("ProFiltTrans_Ran")
data("ProFiltTrans_Shark")
}
\format{
  The format for \code{ProFiltTrans_Ran} is: \cr
List of 7 \cr
 $ Names    : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
 $ maglist  :List of 13 \cr
  ..$ Age0 : num [1:1000, 1:378] 19.9 18.2 20.3 17.9 18.8 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age1 : num [1:1000, 1:378] 22 22.4 22.1 22 21.9 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age2 : num [1:1000, 1:378] 23.1 23.9 23.1 23.6 23.3 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age3 : num [1:1000, 1:378] 24 24.9 24 24.5 24.2 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age4 : num [1:1000, 1:378] 24.9 25.8 24.8 25.5 25.1 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age5 : num [1:1000, 1:378] 25.7 26.6 25.6 26.3 25.9 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age6 : num [1:1000, 1:378] 26.4 27.3 26.2 27.2 26.7 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age7 : num [1:1000, 1:378] 27 27.9 26.9 27.7 27.3 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age8 : num [1:1000, 1:378] 28.1 28.6 28 28.6 28.2 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age9 : num [1:1000, 1:378] 29.2 29.6 29.4 29.8 29.2 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age10: num [1:1000, 1:378] 34 31.7 34.1 32.1 32.7 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age11: num [1:1000, 1:378] 36.5 36.1 36.9 36.2 36.1 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
  ..$ Age12: num [1:1000, 1:378] 40.9 39.6 41.8 39.4 39.9 ... \cr
  .. ..- attr(*, "dimnames")=List of 2 \cr
  .. .. ..$ : chr [1:1000] "result.1" "result.2" "result.3" "result.4" ... \cr
  .. .. ..$ : chr [1:378] "RestUV/Tophat_2700_200.dat" "hst/wfc3/UVIS/f275w.dat" "REST_FRAME/UV2800.dat" "RestUV/Tophat_2800_200.dat" ... \cr
 $ SFHparams: num [1:1000, 1:6] 0.774 0.197 0.978 0.201 0.361 ... \cr
  ..- attr(*, "dimnames")=List of 2 \cr
  .. ..$ : NULL \cr
  .. ..$ : chr [1:6] "m1" "m2" "m3" "m4" ... \cr
 $ cenwave  : num [1:378] 2700 2713 2792 2800 3013 ... \cr
 $ pivwave  : num [1:378] 2699 2708 2791 2800 2992 ... \cr
 $ zsteps   : num [1:14] 0 0.0737 0.1561 0.249 0.3553 ... \cr
 $ ages     : int [1:14] 0 1 2 3 4 5 6 7 8 9 ... \cr
}
\details{
All available filters in \code{ProSpect} that cover u band to the K band are included (373 filters). More specifically we include all filters with a pivot wavelength between 2,500 And and 25,000 Ang. This means almost all useful filter conversions can be computed using this, a sensible age window and \code{\link{filterTranMags}} (see Examples).

The following list items are included:

\describe{
  \item{Names}{Character vector; filter names.}
  \item{maglist}{Numeric list; embedded list of filter magnitudes computed for all filters and lookback times.}
  \item{SFHparams}{Numeric matrix; SFH parameters used for the input galaxies.}
  \item{cenwave}{Numeric vector; central wavelengths for the filters}
  \item{zteps}{Numeric vector; redshifts of the ages samples (nearside of the lookback bin).}
  \item{ages}{Numeric vector; ages of the ages samples in Gyrs (nearside of the lookback bin).}
}
}
\examples{
#Using random SFH samples:

data(ProFiltTrans_Ran)
filterTranMags(ProFiltTrans_Ran$maglist$Age3[,c("g_SDSS", "r_SDSS", "i_SDSS")],
               ProFiltTrans_Ran$maglist$Age3[,"r_VST"])

#Using Shark SFH samples:

data(ProFiltTrans_Shark)
filterTranMags(ProFiltTrans_Shark$maglist$Age3[,c("g_SDSS", "r_SDSS", "i_SDSS")],
               ProFiltTrans_Shark$maglist$Age3[,"r_VST"])
}
\keyword{datasets}
