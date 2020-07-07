\name{plot.ProSpectSED}
\alias{plot.ProSpectSED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ProSpect SED
}
\description{
Simple plot of the various luminosity components in a ProSpectSED object.
}
\usage{
\method{plot}{ProSpectSED}(x, xlim = c(100, 1e+07), ylim = c(100, max(x$StarsUnAtten)),
  xlab = "Wavelength (Ang)", ylab = "Lum (Lsol/Ang)", grid = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class 'ProSpectSED'.
}
  \item{xlim}{
Numeric vector; length 2, x limits.
}
  \item{ylim}{
Numeric vector; length 2, y limits.
}
  \item{xlab}{
Character scalar; x label.
}
  \item{ylab}{
Character scalar; y label.
}
  \item{grid}{
Logical; show background grid?
}
  \item{\dots}{
Passed to \code{magplot}.
}
}
\details{
A simple plot for quick visualisation.
}
\value{
Side effect is to create the plot.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{ProSpectSED}}
}
\examples{
data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}

testSED=ProSpectSED(AGNlum=1e43, tau_birth=1, tau_screen=0.3, tau_AGN=2,
alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, speclib=BC03lr, Dale=Dale_NormTot,
AGN=AGN_UnOb_Sparse, filtout=filtout, Dale_M2L_func=Dale_M2L_func, returnall=TRUE)
plot(testSED)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")

