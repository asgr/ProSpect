\name{SpecLoad}
\alias{SpecLoad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load Generic SSP into ProSpect
}
\description{
Loads generic simple stellar population (SSP) spectral library FITS into ProSpect. This means if users can format a generic spectral library into the required format, then ProSpect can use it once it has been loaded and processed with this function.
}
\usage{
SpecLoad(filename, AgeBins = NULL, Labels = list(Zlab = "Metallicity",
  Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
Character scalar; path to multi-extension spectral library FITS file to be read.
}
  \item{AgeBins}{
Function of numeric vector. If passing a function then the output \option{AgeBins} will be computed as \option{AgeBin}(Age). Alternatively a numeric vector can be passed in, which will be used instead. The resulting length of \option{AgeBin} must be one longer than \option{Age} since it provides the numerical extremes of the bins, where \option{Age} is interpretted as the mid-point. Whether the bins should be treated in linear of log space varies depending on the library (e.g. BC03 is linear, but BPASS is log).
}
  \item{Labels}{
Character vector; vector of plot labels (probably no reason to alter this from the default).
}
}
\details{
Fairly simple function to convert a target external FITS to the internal list format used by \code{\link{ProSpectSED}}.

The extensions and ordering of the target FITS file should look similar to the included BC03lr example used in the Examples below (e.g. note the extension names and the overall ordering): \cr

\tabular{lllll}{
 Ext  \tab       Name  \tab      Class  \tab        Dim  \tab Size.MB \cr
   1  \tab          Z  \tab    numeric  \tab          6  \tab  0.0001 \cr
   2  \tab        Age  \tab    numeric  \tab        221  \tab  0.0017 \cr
   3  \tab    AgeBins  \tab    numeric  \tab        221  \tab  0.0017 \cr
   4  \tab AgeWeights  \tab    numeric  \tab        221  \tab  0.0017 \cr
   5  \tab       Wave  \tab    numeric  \tab       1221  \tab  0.0094 \cr
   6  \tab     Zspec1  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
   7  \tab     Zspec2  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
   8  \tab     Zspec3  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
   9  \tab     Zspec4  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
  10  \tab     Zspec5  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
  11  \tab     Zspec6  \tab     matrix  \tab 221 x 1221  \tab  2.0589 \cr
  12  \tab      Zevo1  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
  13  \tab      Zevo2  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
  14  \tab      Zevo3  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
  15  \tab      Zevo4  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
  16  \tab      Zevo5  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
  17  \tab      Zevo6  \tab data.frame  \tab    221 x 5  \tab  0.0096 \cr
}

The key thing re the multi-extension FITS format is that the \option{Zspec} and \option{Zevo} structures are flattened so there is one matrix (\option{Zpsec}) and one table (\option{Zevo}) per extension. They must have also have extension names starting 'Zspec' and 'Zevo', e.g. 'Zspec1' and 'Zevo1' etc. What follows the leading part of the name is not important. These must be ordered in the same increasing manner as \option{Z}, with \option{Zspec} matrices extensions positioned immediately before \option{Zevo} tables. See the example BC03lr.fits for guidance.
}
\value{
Complex list. See e.g. \code{\link{BC03lr}}. The below is the basic format output:

  \item{Z}{Numeric vector of length Nz describing the availability and location of different metallicity stellar populations.}
  \item{Age}{Numeric vector of length Na giving the mid point of the age in years for the stellar populations.}
  \item{AgeBins}{Numeric vector of length Na+1 giving the age bin extremes in years of the stellar populations. These limits might need to be linearly of log spaced as appropriate.}
  \item{AgeWeights}{Numeric vector of length Na giving the total time in years for the stellar population bins. This is internally calculated as AgeWeights = diff(AgeBins) to conserve mass.}
  \item{Wave}{Numeric vector of length Nw giving the spectral wavelength in Angstroms.}
  \item{Labels}{Character vector of length 5 giving handy axis labels.}
  \item{Zpsec}{List containing 6 (for each metallicity) Na x Nw (i.e. Age rows by Wave columns, so each column should be a contiguous spectrum) matrices that give the full spectral information. The units of the spectra must be Lsol / Ang (for 1 Msol SF). This is pretty standard (i.e. how BC03 comes).}
  \item{Zevo}{List containing 6 (for each metallicity) Na x 5 data frames that contain the evolution information for the stellar material (stellar mass in stars [1] / gas [2] / total [3]; the star formation rate [4] and the stellar mass in remnants [5]). This columns must have the exact names: SMstar / SMgas / SMtot / SFR / SMrem.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{BC03}},\code{\link{EMILES}}, \code{\link{BPASS}}
}
\examples{
\dontrun{
data(BC03lr)
BC03lr_test = SpecLoad(system.file('BC03lr.fits', package="ProSpect"))

#The differences should be very small (parts in 1e-8):
all.equal(BC03lr, BC03lr_test)
}
}
\concept{ speclib }% use one of  RShowDoc("KEYWORDS")
