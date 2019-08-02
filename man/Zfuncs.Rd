\name{Zfunc}
\alias{Zfunc}
\alias{Zfunc_p2}
\alias{Zfunc_massmap_lin}
\alias{Zfunc_massmap_box}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Useful Metallicity Formation Functions
}
\description{
These functions can all be passed into \code{\link{SFHfunc}} and \code{\link{ProSpectSED}}. They should be a good starting point for people interested in constructing their own.
}
\usage{
Zfunc_p2(age, Z1 = 0.02, Z2 = Z1, Z1age = 0, Z2age = Zagemax, Zagemax = 13.8, ...)
Zfunc_massmap_lin(age, Zstart = 1e-4, Zfinal = 0.02, Zagemax = 13.8, massfunc, ...)
Zfunc_massmap_box(age, Zstart = 1e-4, Zfinal = 0.02, yield = 0.03, Zagemax = 13.8,
massfunc, ...)
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
  \item{yield}{
The yield for the closed box model, where Zfinal = Zstart - (yield x ln(gasfrac_final)).
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

\code{Zfunc_massmap_lin} linearly maps the fraction of mass formed for the given SFH against the target start (\option{Zstart}) and finish (\option{Zfinal}) metallicity. Basically this means if you double you fractional mass you will also roughly double the metallicity you form stars with (assuming you start at close to 0 metallicity, which is the default). In detail Zgas = Zstart + (Zfinal - Zstart) x starfrac (where starfrac is the fraction of final stars formed at a given time, i.e. 1 today and 0 at the Big Bang by construction).

\code{Zfunc_massmap_box} maps the fraction of mass formed to Z using Z = Zstart - (yield x ln(gasfrac)). The lowest value for the gas fraction is derived using Zfinal = Zstart - (yield x ln(gasfrac_final)), i.e. gasfrac_final = exp(-(Zfinal - Zstart)/yield). Then Zfinal = Zstart - (yield x ln((1 - starfrac) x gasfrac_final)) (where starfrac is the fraction of final stars formed at a given time, i.e. 1 today and 0 at the Big Bang by construction). Interestingly, if we use this mapping then we naturally predict a gas fraction for the galaxy today!
}
\value{
In all cases the output is the metallicity (Z = 1 - (X + Y), where X is the mass fraction in H and Y is the mass fraction in He) at a given age in years (where 0.02 is solar).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFHfunc}}, \code{\link{ProSpectSED}}, \code{\link{massfunc}}
}
\examples{
curve(Zfunc_p2(age=x, Z1=0.03, Z2=0.01), 0, 13.8e9, xlab='Age/yr', ylab='Z',
ylim=c(0,0.03))
curve(Zfunc_massmap_lin(age=x, Zstart=0.01, Zfinal=0.03, massfunc=massfunc_snorm),
add=TRUE, col='red')
curve(Zfunc_massmap_box(age=x, Zstart=0.01, Zfinal=0.03, massfunc=massfunc_snorm),
add=TRUE, col='blue')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Z }% use one of  RShowDoc("KEYWORDS")
\keyword{ metallicity }% __ONLY ONE__ keyword per line
