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
DEPRECATED - use \code{\link{ProSpectSEDlike}} instead, which provides a more advanced and generic interface for model fitting.

In basic usage, computes the likelihood for a given set of data versus a specified model for the star formation history (SFH), including dust and attentuation. Can also return all of the \code{ProSpect} outputs for a given SFH model produced by \code{SFHp4like}, \code{SFHp5like}, \code{SFHfunclike} and \code{\link{Dale_interp}}.
}
\usage{
SFHp4like(parm = c(8, 9, 10, 10, 0, -0.5, 0.2, -2), Data, massfit = c("burstmass",
"youngmass", "oldmass", "ancientmass"), taufit = c("tau_birth", "tau_screen"),
powfit = NULL, dustfit = c("alpha_SF", "AGNfrac"), zfit = FALSE, verbose = TRUE,
sparse = 5)

SFHp5like(parm = c(8, 9, 10, 10, 10, 0, -0.5, 0.2, -2), Data, massfit = c('burstmass',
'youngmass', 'midmass', 'oldmass', 'ancientmass'), taufit = c('tau_birth', 'tau_screen'),
powfit = NULL, dustfit = c('alpha_SF', 'AGNfrac'), zfit = FALSE, verbose = TRUE,
sparse = 5)

SFHfunclike(parm = c(1, 0, -0.5, 0.2, -2), Data,
massfunc = function(age, SFR = 1){ifelse(age < 1.3e+10, SFR, 0)}, forcemass = FALSE,
unimax = 13.8e9, agescale = 1, massfuncfit = 'SFR', massfuncpos = TRUE,
taufit = c('tau_birth', 'tau_screen'), powfit = NULL, dustfit = c('alpha_SF', 'AGNfrac'),
zfit = FALSE, verbose = TRUE, sparse = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Numeric vector; the current values of the parameters being tested. All standard \code{ProSpect} parameters should be entered as log10, so \option{AGNfrac}=-1 rather than 0.1. The exception is pow_birth and pow_screen, which are in linear units. This encourages better behaviour in all optimisers, since these are all scale not location parameters.
}
  \item{Data}{
Object; a data object containing the observed data to compare the current ProSpect model against. This has to be a list containing a number of vital parts: flux; z; mon.names; parm.names; N; speclib; Dale; filtout; fit; like. Also has optional elements Z; fixed. See Details.
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
Character vector; the names of the Charlot and Fall tau dust attenuation variables to be used (a mixture of tau_birth, tau_screen). See \code{\link{CF_screen}}.
}
  \item{powfit}{
Character vector; the names of the Charlot and Fall power-law dust attenuation variables to be used (a mixture of pow_birth, pow_screen). See \code{\link{CF_screen}}. By default this is not applied, so fitting will not be done to the value of the attenuation power-law.
}
  \item{dustfit}{
Character vector; the names of the dust template variables to be used (a mixture of alpha_SF, AGNfrac).
}
  \item{zfit}{
Logical; whether redshift fitting is being done. If so the first argument of \option{parm} is interpretted to be the redshift, making the stellar mass the second argument if \option{forcemass}=TRUE.
}
  \item{verbose}{
Logical; should verbose output be printed?
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
}
\details{
This function is mostly used for its side effect of producing a likelihood that we can then optimise our model using \code{\link{optim}}. It handles a lot of the awkward data manipulation side which is present when we are using an odd subset of all possible filters etc.

\option{Data} can take most of the options of \code{SFHp4like}, \code{SFHp5like} and \code{SFHfunclike} as inputs to change the defaults, e.g. \option{massfit}, \option{taufit}, \option{powfit}, \option{dustfit}, \option{zfit}, \option{verbose}, \option{sparse}, \option{massfunc}, \option{forcemass}, \option{unimax}, \option{agescale}, \option{massfuncfit}, \option{massfuncpos}.

Some more info on the required \option{Data} list inputs:

\describe{
  \item{Data$flux}{
Data.frame(4 columns, filter name, central wave, flux (Jy), fluxerr (Jy))
}
  \item{Data$prior}{
A prior function that takes in the \option{parm} vector and returns a log-likelihood to be added to LL to form LP. If not provided then it is ignored internally.
}
  \item{Data$speclib}{
Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{Data$Dale}{
Pass in the Dale dust library directly. Must be one of \code{\link{Dale_Orig}}, \code{\link{Dale_Msol}}, \code{\link{Dale_NormTot}}, \code{\link{Dale_NormAGN}}, \code{\link{Dale_NormSFR}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{Data$filtout}{
Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} (see \code{\link{getfilt}}).
}
  \item{Data$fit}{
Character scalar; one of 'optim', 'LD' or 'check'. See Value for full information on what is output in each case.
}
  \item{Data$like}{
Character scalar; specifies the type of likelihood. One of 'norm' (Normal), 'st' (Student-T) or 'chisq' (Chi-Square). If not one of these then all useful output is returned in a list. 
}
  \item{Data$Z}{
Integer vector (optional) of length 4 (\code{SFHp4}), 5 (\code{SFHp5}) or a function (\code{SFHfunc}); specifies which metallicity to use for each star formation phase. The metallicity is the vector location \option{Z} in BC03lr$Z / BC03hr$Z / EMILES$Z., i.e. the default of 5 for BC03lr is close to solar (BC03lr$Z[5]=0.02). For \code{SFHfunc} the functional form takes age as an input (in units of \option{agescale}) and outputs the target metallicity in terms of Fe/H (so 0.02 for solar).
}
  \item{Data$fixed}{
List; optional list of named parameters and values that are being fixed during fitting (i.e. the \code{\link{SFH}} functions will use the fixed values during the fitting process). This includes the redshift z, which will usually be fixed during fitting, i.e. typically Data$fixed$z will need to be specified even if other parameters are not fixed. If this is not specified then vaguely sensible defaults are assumed for the non-fitting arguments (z = 0, burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10, oldmass = 1e+10, ancientmass = 1e+10, tau_birth = 1, tau_screen = 0.3, pow_birth = -0.7, pow_screen = -0.7, alpha_SF = 1.5 and AGNfrac = 0), but the user should really set the fixed values. All fixed values should specified in linear units.
}
}

Outside of the above, any argument that is not \option{parm} and \option{Data} can be passed as part of the named list \option{Data}, rather than having to set it via the individual arguments. This is for two reasons: firstly it means all the user options can be specified in an object that can be saved (to improve reproduceability), but secondly LaplacesDemon can only accept likelihood functions that take exactly two arguments (\option{parm} and \option{Data}), so all options that the user wants to chance have to input as parts of the \option{Data} list structure.
}
\value{
If Data$fit='optim' then the only thing returned is the negative of the log-likelihood as a numeric scalar. Using this output we can minimise our model using \code{\link{optim}}.

If Data$fit='LD' then a list is returned containing LP, Dev, Monitor (log10 dust mass, log10 dust luminosity, log10 AGN luminosity, bolometric dust fraction, bolometric AGN fraction), yhat, parm. Using this output we can optimise our model using \code{LaplacesDemon}.

If Data$fit='check' then lots of things are returned for checking purposes (TBD exactly what). Currently (not exhaustively) includes:

  \item{LP}{Posterior log-likelihood (Data-Model and Prior)}
  \item{LL}{Data-Model log-likelihood}
  \item{Dev}{Deviance: -2.LL}
  \item{FinalSpec}{Numeric matrix; the combined observed flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{FinalPhotom}{Numeric vector; observed photometry in Jansky, as given by \code{\link{Janskycalc}}} in the 'out' column.
  \item{SFH_dust}{List; output of relevant \code{\link{SFH}} function using the current parameter values including dust attenuation.}
  \item{SFH_nodust}{List; output of relevant \code{\link{SFH}} function using the current parameter values but ignoring dust attenuation (effectively \option{tau_birth}=0 and \option{tau_screen}=0.}
  \item{Dust}{Numeric vector; dust properties as determined by \code{\link{dustmass}}.}
  \item{SMtot}{Stellar masses as return by \code{\link{SMstarp4}}, \code{\link{SMstarp5}}, \code{\link{SMstarfunc}}.}
  \item{DaleSpec}{Numeric matrix; the combined observed Dale (dust+AGN) flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{DaleSpec_Dust}{Numeric matrix; the observed Dale dust flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{DaleSpec_AGN}{Numeric matrix; the observed Dale AGN flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{Monitor}{Numeric vector; log10 dust mass, log10 dust luminosity, log10 AGN luminosity, bolometric dust fraction, bolometric AGN fraction.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{SFH}}, \code{\link{magABcalc}}, \code{\link{dustmass}}, \code{\link{optim}}
}

\examples{
## Nothing compact yet, see vignettes.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ likelihood }% use one of  RShowDoc("KEYWORDS")
