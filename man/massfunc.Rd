\name{massfunc}
\alias{massfunc}
\alias{massfunc_const}
\alias{massfunc_p2}
\alias{massfunc_p6}
\alias{massfunc_exp}
\alias{massfunc_CSFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Useful Star Formation Functions
}
\description{
These functions can all be passed into \code{\link{SFHfunc}} and \code{\link{ProSpectSED}}. They should be a good starting point for people interested in constructing their own.
}
\usage{
massfunc_const(age, mSFR = 1, magemax = 13.8e9)
massfunc_p2(age, m1 = 1, m2 = 1, m1age = 0, m2age = magemax, magemax = 13.8e9)
massfunc_p6(age, m1 = 1, m2 = 1, m3 = 1, m4 = 1, m5 = 1, m6 = 1, m1age = 1e5, m2age = 1e8,
  m3age = 1e9, m4age = 5e9, m5age = 9e9, m6age = 1.3e10, magemax = 13.8e9)
massfunc_exp(age, mSFR = 10, mtau = 1, mpivot = magemax, magemax = 13.8e9)
massfunc_CSFH(age, mSFR = 10, mpeak = 1e10, mperiod = 1e9, mskew = 0.5, magemax = 13.8e9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
Numeric vector; age in lookback time (yr).
}
  \item{mSFR}{
Numeric scalar; SFR scaling (Msol/yr). For most functions this ends up being the peak SFR [\code{massfunc_const} / \code{massfunc_exp} / \code{massfunc_CSFH}].  
}
  \item{m1}{
Numeric scalar; SFR at 1st node (Msol/yr) [\code{massfunc_p2} / \code{massfunc_p6}].
}
  \item{m2}{
Numeric scalar; SFR at 2nd node (Msol/yr) [\code{massfunc_p2} / \code{massfunc_p6}].
}
  \item{m3}{
Numeric scalar; SFR at 3rd node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m4}{
Numeric scalar; SFR at 4th node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m5}{
Numeric scalar; SFR at 5th node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m6}{
Numeric scalar; SFR at 6th node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m1age}{
Numeric scalar; age at 1st node (yr) [\code{massfunc_p2} / \code{massfunc_p6}].
}
  \item{m2age}{
Numeric scalar; age at 2nd node (yr) [\code{massfunc_p2} / \code{massfunc_p6}].
}
  \item{m3age}{
Numeric scalar; age at 3rd node (yr) [\code{massfunc_p6}].
}
  \item{m4age}{
Numeric scalar; age at 4th node (yr) [\code{massfunc_p6}].
}
  \item{m5age}{
Numeric scalar; age at 5th node (yr) [\code{massfunc_p6}].
}
  \item{m6age}{
Numeric scalar; age at 6th node (yr) [\code{massfunc_p6}].
}
  \item{mtau}{
Numeric scalar; exponential decline parameter tau [\code{massfunc_exp}].
}
  \item{mpivot}{
Numeric scalar; exponential pivot age where SFR passes though \option{mSFR} (yr) [\code{massfunc_exp}].
}
  \item{mpeak}{
Numeric scalar; location of peak SFR (yr) [\code{massfunc_CSFH}].
}
  \item{mperiod}{
Numeric scalar; time period of the star formation, the standard deviation of the Normal when \option{mskew}=0 (yr) [\code{massfunc_CSFH}].
}
  \item{mskew}{
Numeric scalar; the skew of the star formation history, 0 means perfectly Normal, -ve means long tail to old ages, +ve means long tail to young ages [\code{massfunc_CSFH}].
}
  \item{magemax}{
Numeric scalar; max age allowed, SFR is set to 0 beyond this value (yr) [all].
}
}
\details{
Below are the brief functional forms.

\code{massfunc_const}: 

SFR = ifelse(age < magemax, mSFR, 0)

\code{massfunc_exp}:

SFR = mSFR*exp(-mtau*((mpivot-age)/mpivot)) [also 0 outside magemax]

\code{massfunc_p2}:

SFR = m1 + ((m2-m1)/(m2age-m1age))*(age[sel]-m1age) [also 0 outside of m1age / m2age / magemax]

\code{massfunc_p6}:

SFR = splinefun(log10(c(m1age,m2age,m3age,m4age,m5age,m6age)), c(m1,m2,m3,m4,m5,m6), method='monoH.FC')(log10(age) [also 0 outside magemax]

\code{massfunc_CSFH}:

SFR = mSFR*dnorm(((age-mpeak)/mperiod)*(exp(mskew))^asinh((age-mpeak)/mperiod))*sqrt(2*pi) [also 0 outside magemax]
}
\value{
In all cases the output is the star formation rate (SFR) at a given age in years.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFHfunc}}, \code{\link{ProSpectSED}}
}
\examples{
curve(massfunc_const, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p2, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p6, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_exp, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_CSFH, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SFR }% use one of  RShowDoc("KEYWORDS")
