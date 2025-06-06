\name{SKIRTOR_interp}
\alias{SKIRTOR_interp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate the SKIRTOR AGN Templates
}
\description{
Interpolates between non-computed values of the Stalevski SKIRTOR AGN Templates using a mixture of schemes. This is similar to the approach used for \code{\link{Dale_interp}}.
}
\usage{
SKIRTOR_interp(lum = 1e+44, ta = 1, p = 1, q = 1, ct = 40, rm = 60,
                             an=30, SKIRTOR = NULL)
}
\arguments{
  \item{lum}{
Numeric scalar; AGN bolometric luminosity in erg/s. This is of the illuminating accretion source, so actual observed output will vary depending on the geometry of the dust torus.
}
  \item{ta}{
Numeric scalar; Average edge-on optical depth at 9.7 micron.
}
  \item{p}{
Numeric scalar; Index for radial density gradient.
}
  \item{q}{
Numeric scalar; Index for angular density gradient.
}
  \item{ct}{
Numeric scalar; Opening angle of torus in degrees.
}
  \item{rm}{
Numeric scalar; Ratio of outer radius to inner radius.
}
  \item{an}{
Numeric scalar; angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned (this is opposite to the definition of the SKIRTOR model itself, in order to align with the Fritz model).
}
  \item{SKIRTOR}{
The \option{SKIRTOR} data. Providing this directly reduces disk IO since it will use the copy already loaded into RAM. If not provided it will load the data from disk.
}

}
\details{
The below is taken directly from the included documentation:

We modeled the distribution of the dust in the torus as a two-phase medium, consisting of a large number of high-density clumps embedded in a smooth dusty component of low density. 
Previously, we have found that an advantageous property of such two-phase medium is that it can produce attenuated silicate features, while at the same time have a pronounced NIR emission, which is challenging for both smooth and clumpy models alone.


}
\value{
Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed luminosity (erg/s/Ang).
}
\references{
Stalevski et al, 2012, MNRAS, 420, 2756 \cr
Stalevski et al, 2016, MNRAS, 458, 2288
}
\author{
Adam Marshall
}
\seealso{
\code{\link{SKIRTOR}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SKIRTOR }% use one of  RShowDoc("KEYWORDS")
