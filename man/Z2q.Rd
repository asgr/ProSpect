\name{Z2q}
\alias{Z2q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metal to Ionisation Conversion
}
\description{
Uses the Z[gas] to ionisation conversion function suggested in Orsi 2014/
}
\usage{
Z2q(Zgas = 0.02, q0 = 2.8e+07, g0 = -1.3, Z0 = 0.012)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Zgas}{
Numeric scaler; gas phase metallicity.
}
  \item{q0}{
Numeric scaler; reference ionisation parameter. The ionization parameter is the ratio of the incident ionizing photon flux (S[H0]) to the gas density (n), q=S[H0]/n.
}
  \item{g0}{
Numeric scalar; power law scaling of relationship.
}
  \item{Z0}{
Numeric scaler; reference Z[gas], where gas with metallicity \option{Z0} has ionisation parameter \option{q0}.
}
}
\details{
Mapping given by: q(Z) = q0(Z/Z0)^g0.

This function probably will not be used by people directly, but it is used inside of \code{\link{emissionLines}} to make a sensible mapping when \option{q} is not specified.
}
\value{
Output is numeric scalar value of \option{q} with units cm/s.
}
\references{
Orsi et al, 2014, MNRAS, 443, 799
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{emissionLines}}
}
\examples{
# For solar metallicity:

Z2q(0.02)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ metallicity }% use one of  RShowDoc("KEYWORDS")
