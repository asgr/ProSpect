\name{massfunc}
\alias{massfunc}
\alias{massfunc_const}
\alias{massfunc_p2}
\alias{massfunc_p3}
\alias{massfunc_p3_burst}
\alias{massfunc_p4}
\alias{massfunc_p6}
\alias{massfunc_b5}
\alias{massfunc_exp}
\alias{massfunc_exp_burst}
\alias{massfunc_CSFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Useful Star Formation Functions
}
\description{
These functions can all be passed into \code{\link{SFHfunc}} and \code{\link{ProSpectSED}}. They should be a good starting point for people interested in constructing their own.
}
\usage{
massfunc_const(age, mSFR = 1, magemax = 13.8)
massfunc_p2(age, m1 = 1, m2 = m1, m1age = 0, m2age = magemax, magemax = 13.8)
massfunc_p3(age, m1 = 1, m2 = m1, m3 = m2, m1age = 1e-4, m2age = 7, m3age = 13,
  magemax = 13.8)
massfunc_p3_burst(age, mburst = 0, m1 = 1, m2 = m1, m3 = m2, mburstage = 0.1,
  m1age = 1e-4, m2age = 7, m3age = 13, magemax = 13.8)
massfunc_p4(age, m1 = 1, m2 = m1, m3 = m2, m4 = m3, m1age = 1e-4, m2age = 2,
  m3age = 9, m4age = 13, magemax = 13.8)
massfunc_p6(age, m1 = 1, m2 = m1, m3 = m2, m4 = m3, m5 = m4, m6 = m5, m1age = 1e-4,
  m2age = 0.1, m3age = 1, m4age = 5, m5age = 9, m6age = 13, magemax = 13.8)
massfunc_b5(age, m1 = 1, m2 = m1, m3 = m3, m4 = m3, m5 = m4, m1age = 0, m2age = 0.1,
  m3age = 1, m4age = 5, m5age = 9, m6age = 13, magemax = 13.8)
massfunc_exp(age, mSFR = 10, mtau = 1, mpivot = magemax, magemax = 13.8)
massfunc_exp_burst(age, mburst = 0, mSFR = 10, mburstage = 0.1, mtau = 1,
  mpivot = magemax, magemax = 13.8)
massfunc_CSFH(age, mSFR = 10, mpeak = 10, mperiod = 1, mskew = 0.5, magemax = 13.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
Numeric vector; age in lookback time where 0 means now (yr).
}
  \item{mSFR}{
Numeric scalar; SFR scaling (Msol/yr). For most functions this ends up being the peak SFR [\code{massfunc_const} / \code{massfunc_exp} / \code{massfunc_exp_burst} / \code{massfunc_CSFH}].  
}
  \item{mburst}{
Numeric scalar; SFR of a recent Age < \option{mburstage} yr burst (Msol/yr) [\code{massfunc_p3_burst} / \code{massfunc_exp_burst}].
}
  \item{mburstage}{
Numeric scalar; max age of recent SFR burst [\code{massfunc_p3_burst} / \code{massfunc_exp_burst}].
}
  \item{m1}{
Numeric scalar; SFR at 1st node (Msol/yr) [\code{massfunc_p2} / \code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p4} / \code{massfunc_p6}].
}
  \item{m2}{
Numeric scalar; SFR at 2nd node (Msol/yr) [\code{massfunc_p2} / \code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p4} / \code{massfunc_p6}].
}
  \item{m3}{
Numeric scalar; SFR at 3rd node (Msol/yr) [\code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p6} / \code{massfunc_p4}].
}
  \item{m4}{
Numeric scalar; SFR at 4th node (Msol/yr) [\code{massfunc_p6} / \code{massfunc_p4}].
}
  \item{m5}{
Numeric scalar; SFR at 5th node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m6}{
Numeric scalar; SFR at 6th node (Msol/yr) [\code{massfunc_p6}].
}
  \item{m1age}{
Numeric scalar; age at 1st node (Gyr) [\code{massfunc_p2} / \code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p4} / \code{massfunc_p6}].
}
  \item{m2age}{
Numeric scalar; age at 2nd node (Gyr) [\code{massfunc_p2} / \code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p4} / \code{massfunc_p6}].
}
  \item{m3age}{
Numeric scalar; age at 3rd node (Gyr) [ \code{massfunc_p3} / \code{massfunc_p3_burst} / \code{massfunc_p6} / \code{massfunc_p4}].
}
  \item{m4age}{
Numeric scalar; age at 4th node (Gyr) [\code{massfunc_p6} / \code{massfunc_p4}].
}
  \item{m5age}{
Numeric scalar; age at 5th node (Gyr) [\code{massfunc_p6}].
}
  \item{m6age}{
Numeric scalar; age at 6th node (Gyr) [\code{massfunc_p6}].
}
  \item{mtau}{
Numeric scalar; exponential decline parameter tau [\code{massfunc_exp} / \code{massfunc_exp_burst}].
}
  \item{mpivot}{
Numeric scalar; exponential pivot age where SFR passes though \option{mSFR} (Gyr) [\code{massfunc_exp} / \code{massfunc_exp_burst}].
}
  \item{mpeak}{
Numeric scalar; location of peak SFR (Gyr) [\code{massfunc_CSFH}].
}
  \item{mperiod}{
Numeric scalar; time period of the star formation, the standard deviation of the Normal when \option{mskew}=0 (Gyr) [\code{massfunc_CSFH}].
}
  \item{mskew}{
Numeric scalar; the skew of the star formation history, 0 means perfectly Normal, -ve means long tail to old ages, +ve means long tail to young ages [\code{massfunc_CSFH}].
}
  \item{magemax}{
Numeric scalar; max age allowed, SFR is set to 0 beyond this value (Gyr) [all].
}
}
\details{
Below are the brief functional forms.

\code{massfunc_const} produces constant star formation until a certain age: 

SFR = ifelse(age < magemax, mSFR, 0)

\code{massfunc_exp} produces exponentially declining star formation from a certain age:

SFR = mSFR*exp(-mtau*((mpivot-age)/mpivot)) [also 0 above magemax]

\code{massfunc_p2} produces linearly varying star formation between 2 age nodes:

SFR = m1 + ((m2-m1)/(m2age-m1age))*(age-m1age) [also 0 outside of m1age / m2age / magemax]

\code{massfunc_p3} produces a smooth spline star formation that passes through 3 age nodes:

SFR = splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age) [also 0 above magemax]

\code{massfunc_p3_burst} produces a smooth spline star formation that passes through 3 age nodes and has a recent burst:

SFR = splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age) [also 0 above magemax]
SFR[age<mburstage] = SFR[age<mburstage] + mburst

\code{massfunc_p4} produces a smooth spline star formation that passes through 6 age nodes:

SFR = splinefun(log10(c(m1age,m2age,m3age,m4age)), c(m1,m2,m3,m4), method='monoH.FC')(log10(age) [also 0 above magemax]

\code{massfunc_p6} produces a smooth spline star formation that passes through 6 age nodes:

SFR = splinefun(log10(c(m1age,m2age,m3age,m4age,m5age,m6age)), c(m1,m2,m3,m4,m5,m6), method='monoH.FC')(log10(age) [also 0 above magemax]

\code{massfunc_CSFH} produces a rising and falling star formation that looks a bit like the CSFH (can be perfectly Normal if \option{mskew}=0):

SFR = mSFR*dnorm(((age-mpeak)/mperiod)*(exp(mskew))^asinh((age-mpeak)/mperiod))*sqrt(2*pi) [also 0 above magemax]
}
\value{
In all cases the output is the star formation rate (Msol/yr) at a given age in years.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{SFHfunc}}, \code{\link{ProSpectSED}}, \code{\link{Zfunc}}
}
\examples{
curve(massfunc_const, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p2, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p3, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p3_burst, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p4, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_p6, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_exp, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
curve(massfunc_CSFH, 0, 14e9, xlab='Age/yr', ylab='SFR / Msol/yr')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SFR }% use one of  RShowDoc("KEYWORDS")