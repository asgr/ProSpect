\name{SFHlike}
\alias{SFHlike}
\alias{SFHp4like}
\alias{SFHp5like}
\alias{SFHfunclike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Star Formation History Likelihood
}
\description{
In basic usage, computes the likelihood for a given set of data versus a specified model for the star formation history (SFH), including dust and attentuation. Can also return all of the \code{ProSpect} outputs for a given SFH model produced by \code{\link{SFHp4}} and \code{\link{Dale_interp}}.
}
\usage{
SFHp4like(parm = c(8, 9, 10, 10, 1, 0.3, 1.5, 0), Data, massfit = c("burstmass",
"youngmass", "oldmass", "ancientmass"), taufit = c("tau_birth", "tau_screen"),
dustfit = c("alpha_SF", "AGNfrac"), zfit = FALSE, filters = c("FUV", "NUV", "u_SDSS", "g_SDSS",
"r_SDSS", "i_SDSS", "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2",
"W3", "W4", "P100", "P160", "S250", "S350", "S500"), verbose = TRUE, like = TRUE,
speclib = NULL, Dale = NULL, filtout = NULL, sparse = 1)

SFHp5like(parm = c(8, 9, 10, 10, 10, 1, 0.3, 1.5, 0), Data, massfit = c('burstmass',
'youngmass', 'midmass', 'oldmass', 'ancientmass'), taufit = c('tau_birth', 'tau_screen'),
dustfit = c('alpha_SF', 'AGNfrac'), zfit = FALSE, filters = c('FUV', 'NUV', 'u_SDSS',
'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1',
'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose = TRUE,
like = TRUE, speclib = NULL, Dale = NULL, filtout = NULL, sparse = 1)

SFHfunclike(parm = c(1, 1, 0.3, 1.5, 0), Data, massfunc = function(age, SFR = 1)
{ ifelse(age < 1e+10, SFR, 0)}, forcemass = FALSE, unimax = 13.8e9, agescale = 1,
massfuncfit = c('SFR'), massfuncpos = TRUE, taufit = c('tau_birth', 'tau_screen'),
dustfit = c('alpha_SF', 'AGNfrac'), zfit = FALSE, filters=c('FUV', 'NUV', 'u_SDSS',
'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1',
'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose=TRUE, like=TRUE,
speclib=NULL, Dale=NULL, filtout=NULL, sparse = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Numeric vector; the current values of the parameters being tested.
}
  \item{Data}{
Object; a data object containing the observed data to compare the current ProSpect model against.
}
  \item{massfunc}{
Function; a function which will take \code{age} (in Yrs) as an input and return SFR (in Msun/Yr). Can be passed further arguments via \dots (this is useful if you do not want to hard code values into a function, i.e. for fitting purposes).
}
  \item{forcemass}{
Logical; if FALSE then the stellar mass is directly computed from the \option{massfunc} function, i.e. the output of this must be a true star formation rate in Msun per year. If TRUE then the first argument of \option{parm} (or second if \option{zfit}=TRUE) is interpretted to be the target stellar mass in log10 units, and the output of \option{massfunc} will linearly scaled as required to achieve this.
}
  \item{unimax}{
Numeric scalar; maximum allowed age of any stellar population relative to z=0 (i.e. today). Any star formation implied by \option{massfunc} that is older than this will be set to 0.
}
  \item{agescale}{
Numeric scalar; scaling to apply to the ages (which will passed in years) to \option{massfunc}. E.g. if left at \option{agescale}=1 then \option{massfunc} is expecting the age in years, but if set to 1e-6/1e-9 then \option{massfunc} is expecting the age in Myrs/Gyrs.
}
  \item{massfit}{
Character vector; the names of the mass variables to be used (a mixture of burstmass, youngmass, oldmass, ancientmass).
}
  \item{massfuncfit}{
Character vector; the names of the mass function parameters in \option{massfunc} that exist at the front of \option{parm}.
}
\item{massfuncpos}{
Logical vector; flags whether certain \option{massfuncfit} parameters provided by \option{parm} should be positive only in \option{massfunc}. This is a simple way to stop certain functions from failing during sampling where it might go negative incorrectly (e.g. standard deviation). For more sophisticated behaviour the keen user will need to write their own likelihood function.
}
  \item{taufit}{
Character vector; the names of the tau dust attenuation variables to be used (a mixture of tau_birth, tau_screen).
}
  \item{dustfit}{
Character vector; the names of the dust template variables to be used (a mixture of alpha_SF, AGNfrac).
}
  \item{zfit}{c(21)
Logical; whether redshift fitting is being done. If so the first argument of \option{parm} is interpretted to be the redshift, making the stellar mass the second argument if \option{forcemass}=TRUE.
}
  \item{filters}{
Character vector; names fo the filters to use for comparison.
}
  \item{verbose}{
Logical; should verbose output be printed?
}
  \item{like}{
Logical; should only the likelihood be output (TRUE, useful when optimising), or should all useful output be produced in a list (FALSE). See Details for full information on what is output in each case.
}
  \item{speclib}{
Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{Dale}{
Pass in the Dale dust library directly. Must be one of \code{\link{Dale_Orig}}, \code{\link{Dale_Msol}}, \code{\link{Dale_NormTot}}, \code{\link{Dale_NormAGN}}, \code{\link{Dale_NormSFR}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{filtout}{
Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} (see \code{\link{getfilt}}).
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
}
\details{
This function is mostly used for its side effect of producing a likelihood that we can then optimise our model using \code{\link{optim}}. It handles a lot of the awkward data manipulation side which is present when we are using an odd subset of all possible filters etc.
}
\value{
If \option{like}=TRUE (the default) then the only thing returned is the log-likelihood as a numeric scalar. Using this output we can optimise our model using \code{\link{optim}}.

If \option{like}=FALSE then lots of things are returned (TBD exactly what). Currently includes:

  \item{like}{log-likelihood}
  \item{FinalSpec}{Numeric matrix; the combined observed flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{FinalPhotom}{Numeric vector; observed photometry in Jansky, as given by \code{\link{Janskycalc}}} in the \code{out} column.
  \item{SFH_dust}{List; output of \code{\link{SFHp4}} using the current parameter values including dust attenuation.}
  \item{SFH_nodust}{List; output of \code{\link{SFHp4}} using the current parameter values but ignoring dust attenuation (effectively \option{tau_birth}=0 and \option{tau_screen}=0.}
  \item{Dust}{Numeric vector; dust properties as determined by \code{\link{dustmass}}.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{SFHp4}}, \code{\link{magABcalc}}, \code{\link{dustmass}}, \code{\link{optim}}
}

\examples{
## Nothing compact yet, see vignettes.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }% use one of  RShowDoc("KEYWORDS")
