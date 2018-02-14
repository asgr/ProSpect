\name{SFHlikefunc}
\alias{SFHlikefunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Star Formation History Likelihood
}
\description{
In basic usage, computes the likelihood for a given set of data versus a specified model for the star formation history (SFH), including dust and attentuation. Can also return all of the \code{ProSpect} outputs for a given SFH model produced by \code{\link{SFHp4}} and \code{\link{Dale_interp}}.
}
\usage{
SFHlikefunc(parm = c(8, 9, 9, 9, 1, 0.3, 1.5, 0), Data, massfit = c("burstmass",
"youngmass", "oldmass", "ancientmass"), taufit = c("tau_birth", "tau_screen"),
dustfit = c("alpha_SF", "AGNfrac"), filters = c("FUV", "NUV", "u_SDSS", "g_SDSS",
"r_SDSS", "i_SDSS", "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "W1", "W2",
"W3", "W4", "P100", "P160", "S250", "S350", "S500"), verbose = TRUE, like = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
Numeric vector; the current values of the parameters being tested.
}
  \item{Data}{
Object; a data object containing the observed data to compare the current ProSpect model against.
}
  \item{massfit}{
Character vector; the names of the mass variables to be used (a mixture of burstmass, youngmass, oldmass, ancientmass).
}
  \item{taufit}{
Character vector; the names of the tau dust attenuation variables to be used (a mixture of tau_birth, tau_screen).
}
  \item{dustfit}{
Character vector; the names of the dust template variables to be used (a mixture of alpha_SF, AGNfrac).
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
}
\details{
This function is mostly used for its side effect of producing a likelihood that we can then optimise our model using \code{\link{optim}}. It handles a lot of the awkward data manipulation side which is present when we are using an odd subset of all possible filters etc.
}
\value{
If \option{like}=TRUE (the default) then the only thing returned is the log-likelihood as a numeric scalar. Using this output we can optimise our model using \code{\link{optim}}.

If \option{like}=FALSE then lots of things are returned (TBD exactly what). Currently includes:

  \item{like}{log-likelihood}
  \item{FinalSpec}{Numeric matrix; the combined observed flux spectrum (column 1 is wavelength; column 2 is the observed flux in erg/s / cm^2 / Ang.}
  \item{FinalPhotom}{Numeric matrix; observed photometry in AB mag, as given by \code{\link{magABcalc}}}.
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
