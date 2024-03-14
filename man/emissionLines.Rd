\name{emissionLines}
\alias{emissionLines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Emission Line Generation
}
\description{
Code to generate emission line specra based on Levesque, Kewley & Larson (2010).
}
\usage{
emissionLines(Ha_lum=NULL, Hb_lum=NULL, Hlines_lum=NULL, All_lum=NULL, SFR=NULL,
  Z = 0.02, q = NULL, veldisp = 50, LSF = NULL, z_LSF = 0, lumscale = 21612724,
  log = TRUE, range = 5, res = 0.5, LKL10 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ha_lum}{
Numeric scalar; scaling to be applied to Halpha line, all other lines are scaled relative to this. Uses/requires \code{\link{LKL10_NormHalpha}} if \option{LKL10}=NULL.
}
  \item{Hb_lum}{
Numeric scalar; scaling to be applied to Hbeta line, all other lines are scaled relative to this. Uses/requires \code{\link{LKL10_NormHbeta}} if \option{LKL10}=NULL.
}
  \item{Hlines_lum}{
Numeric scalar; scaling to be applied all Hydrogen lines, all other lines are scaled relative to this. Uses/requires \code{\link{LKL10_NormHlines}} if \option{LKL10}=NULL.
}
  \item{All_lum}{
Numeric scalar; scaling to be applied all emission lines. Uses/requires \code{\link{LKL10_NormAll}} if \option{LKL10}=NULL.
}
  \item{SFR}{
Numeric scalar; scaling to be applied to Halpha line using \code{\link{SFR2Lum}} function to scale via the star formation rate, all other lines are scaled relative to this. Uses/requires \code{\link{LKL10_NormHalpha}} if \option{LKL10}=NULL.
}
  \item{Z}{
Numeric scalar; metallicity. If \option{log}=TRUE this will be interpolated in log space if it falls between templates.
}
  \item{q}{
Numeric scalar; radiation field. If \option{log}=TRUE this will be interpolated in log space if it falls between templates.
}
  \item{veldisp}{
Numeric scalar; velocity dispersion (km/s).
}
  \item{LSF}{
Function, numeric matrix/data.frame, numeric scalar; describes the instrumental line spread function is km/s. This should either scalar (constant throughout), a two column matrix/data.frame (wavelength in Ang, LSF in km/s) or a function that take wavelength (Ang) as an input and returns LSF in km/s. NULL means this is ignored (effectively same as setting to 0, but more efficient).
}
  \item{z_LSF}{
Numeric scalar; the redshift the LSF is effectively provided at. This is because instruments provide and observed frame LSF, so to compute what it is in our intrinsic frame we need to know the redshift we are working at. If \option{LSF} is NULL, this is not required.  
}
  \item{lumscale}{
Numeric scalar; luminosity scaling to use in \code{\link{SFR2Lum}}. Unless you know what you are doing, leave at the default.
}
  \item{log}{
Logical; determines if the interpolation is done in log space when using \code{\link{interp_param}}.
}
  \item{range}{
Numeric scalar; how many sigma either side of the line should the velocity dispersion spread flux be evaluated (starts to become inaccurate if lower than 5).
}
  \item{res}{
Numeric scalar; the relative resolution of the line evaluation in fractions of \option{sigma}.
}
  \item{LKL10}{
List; LKL template to use. If \option{LKL10}=NULL then the appropriate template is loaded internally (this is the safest option). To speed things up the user can pass in the appropriate \code{\link{LKL10}} template, but care must be taken to ensure the scaling used makes use of the correct dataset.
}
}
\details{
Can be used by users directly (perhaps to generate your own line features to add onto a target spectrum), or more generally.
}
\value{
Data.frame with two columns:
\item{wave }{Wavelength (Angstroms).}
\item{lum }{Emission features scaled as requested (Lum/Ang).}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{LKL10}}, \code{\link{SFR2Lum}}, \code{\link{Z2q}}, \code{\link{interp_param}}
}
\references{
Levesque, Kewley & Larson, 2010, AJ, 139, 712
}
\examples{
plot(emissionLines(Ha_lum=1, veldisp=100), type='l', xlim=c(6400,6800), ylim=c(0,0.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ emission }% use one of  RShowDoc("KEYWORDS")
