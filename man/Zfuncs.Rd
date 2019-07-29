\name{Zfunc}
\alias{Zfunc}
\alias{Zfunc_p2}
\alias{Zfunc_massmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Useful Metallicity Formation Functions
}
\description{
These functions can all be passed into \code{\link{SFHfunc}} and \code{\link{ProSpectSED}}. They should be a good starting point for people interested in constructing their own.
}
\usage{
Zfunc_p2(age, Z1 = 0.02, Z2 = Z1, Z1age = 0, Z2age = Zagemax, Zagemax = 13.8, ...)
Zfunc_massmap(age, Zstart = 1e-4, Zfinal = 0.02, Zagemax = 13.8, massfunc, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
Numeric vector; age in lookback time where 0 means now (yr).
}
  \item{Z1}{
Numeric scalar; Z at 1st node (Msol/yr) [\code{Zfunc_p2}].
}
  \item{Z2}{
Numeric scalar; Z at 2nd node (Msol/yr) [\code{Zfunc_p2}].
}
  \item{Z1age}{
Numeric scalar; age at 1st node (Gyr) [\code{Zfunc_p2}].
}
  \item{Z2age}{
Numeric scalar; age at 2nd node (Gyr) [\code{Zfunc_p2}].
}
  \item{Zagemax}{
Numeric scalar; max age allowed, Z is set to 1e-04 beyond this value (Gyr) [all].
}
  \item{Zstart}{
Numeric scalar; Z at the start of star formation (Msol/yr) [\code{Zfunc_massmap}].
}
  \item{Zfinal}{
Numeric scalar; Z at the end of star formation (Msol/yr) [\code{Zfunc_massmap}].
}
  \item{massfunc}{
Function; mass function to have metallicity mapped back against [\code{Zfunc_massmap}], see \code{\link{massfunc}}.
}
  \item{\dots}{
All further arguments are passed directly into \option{massfunc}, so by this mechanism you can set function specific arguments that control the star formation history etc [\code{Zfunc_massmap}].
}
}
\details{
Below are the brief functional forms.

\code{Zfunc_p2} produces linearly varying star formation between 2 age nodes:

Z = Z1 + ((Z2-Z1)/(Z2age-Z1age))*(age-Z1age) [also 1e-4 outside of Z1age / Z2age / Zagemax]

\code{Zfunc_massmap} linearly maps the fraction of mass formed for the given SFH against the target start (\option{Zstart}) and finish (\option{Zfinal}) metallicity. Basically this means if you double you fractional mass you will also roughly double the metallicity you form stars with (assuming you start at close to 0 metallicity, which is the default).
}
\value{
In all cases the output is the metallicity (Z=Fe/H) at a given age in years (where 0.02 is solar).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFHfunc}}, \code{\link{ProSpectSED}}, \code{\link{massfunc}}
}
\examples{
curve(Zfunc_p2, 0, 14e9, xlab='Age/yr', ylab='Z')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Z }% use one of  RShowDoc("KEYWORDS")
\keyword{ metallicity }% __ONLY ONE__ keyword per line
