\name{speclib_download}
\alias{speclib_download}
\alias{speclib_FITSload}
\alias{speclib_check}
\alias{speclib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProSpect Stellar Population Library Downloader
}
\description{
Simple interface to show all the currently available ProSpect SSPs, and download them as needed.
}
\usage{
speclib_download(stellpop = "avail", destpath = "", type = 'rda',
  URL = "https://tinyurl.com/prospect-speclib/")

speclib_FITSload(file, Labels = list(Zlab = "Metallicity",
  Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"),
  check = FALSE)

speclib_check(speclib, structure = TRUE, coverage = TRUE, Labels = list(
  Zlab = "Metallicity", Agelab = "Time since ZAM / Yrs", Wavelab = "Wavelength / Ang",
  Lumlab = "Lsun / Ang (for 1 Msun SF)", LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stellpop}{
Character scalar; either 'avail' which will then display all the SSPs currently available at \option{URL} and other key information, or the name of the target SSP (that must be available!). Note BC03lr_Pd_Ch / BC03hr_Pd_Ch / EMILES_Bs_Ch / BPASS_St_Ch are in the core \code{ProSpectData} package, but they are also included here for completeness (and in case in future some are removed from the core package). The equivalent package version are named BC03lr, BC03hr, EMILES and BPASS respectively.
}
  \item{type}{
Character scalar; type of spectral library to download. Either 'rda' (default) or 'fits'. Note in future perhaps only one of the rda/fits SSP type will be available for a given SSP. This should be clear from checking the \option{stellpop} = 'avail' information.
}
  \item{destpath}{
Character scalar; where the target ProSpect SSP should be downloaded to.
}
  \item{URL}{
Character scalar; location of the ProSpect SSP libraries. This should not need to be altered by the user.
}
  \item{file}{
Character scalar; full path and name of target FITS spectral library to load.
}
  \item{Labels}{
List; list of plot labels (probably no reason to alter this from the default). These will either be created as written, or checked.
}
  \item{check}{
Logical; should checks be run on the spectral library being returned? Uses \code{speclib_check}. If checks fail the SSP is still returned, and can then be checked and potentially fixed by the user.
}
  \item{speclib}{
List; spectral library to pass through \code{ProSpect} format checks.
}
  \item{structure}{
Logical; should the structure of \option{speclib} be checked. If TRUE, results are printed to screen.
}
  \item{coverage}{
Logical; should the coverage (Z, Age, Wave) of \option{speclib} be checked? If TRUE, results are printed to screen.
}
}
\details{
This interface has proven necessary in order to provide flexible SSPs without making the core package too big. We have also had issue with the package getting close to the file size limit of GitHub. In future all additional SSP will be provided through this mechanism.

Internally a \code{ProSpect} compatible format (e.g. as seen for \code{\link{BC03}} / \code{\link{EMILES}} / \code{\link{BPASS}}, as output by \code{speclib_FITSload}, and as verified with \code{speclib_check}) looks like a list containing:

Z: Numeric vector of available template metallicities (ascending order).
Age: Numeric vector of available template ages in Yrs (ascending order).
AgeBins: Numeric vector of available template age bin limits in Yrs (ascending order, must be 1 longer than "Age" vector to cover lower and upper limits).
AgeWeights: Numeric vector of available template age weights (bin sizes) in Yrs (same length as "Age").
Wave: Numeric vector of spectral wavelength in Angstroms (ascending order).
Labels: List of the labels required for plotting (see \code{\link{BC03}}).
Zspec: List containing the spectra in units Lsun / Ang (for 1 Msun SF). The list must correspond to the 'Z' vector metallicities. Each list item is a numeric matrix of spectra with Nrow equal to the number of age bins (corresponding with 'Age') and Ncol the number of wavelength element (corresponding with 'Wave').
Zevo: List containing the mean evolutionary tracks of the SSPs (as determined by the route isochrones). The list must correspond to the 'Z' vector metallicities. Each list item is a data.frame with columns 'SMstar' (fraction of mass in luminous stars) 'SMgas' (fraction of mass in gas) 'SMtot' (fraction of total mass, always 1 in practice) 'SFR' (star formation rate, always 0 in practice) 'SMrem' (fraction of mass in stellar remnants, e.g. black holes etc). Nrow should be equal to the number of age bins (corresponding with 'Age').

For FITS loading, the required format is a multi-extension FITS file containing (and named via EXTNAME):

Z: Numeric vector of available template metallicities (ascending order). The length of this vector is N_Z (as referred to below).
Age: Numeric vector of available template ages in Yrs (ascending order).
AgeBins: Numeric vector of available template age bin limits in Yrs (ascending order, must be 1 longer than "Age" vector to cover lower and upper limits).
AgeWeights: Numeric vector of available template age weights (bin sizes) in Yrs (same length as "Age").
Wave: Numeric vector of spectral wavelength in Angstroms (ascending order).
Zspec[1:N_Z]: Extensions containing the spectra in units Lsun / Ang (for 1 Msun SF). Each numbered extension must correspond to the 'Z' vector metallicities. Each extension is a numeric matrix of spectra with Nrow equal to the number of age bins (corresponding with 'Age') and Ncol the number of wavelength element (corresponding with 'Wave').
Zevo[1:N_Z]: Extensions containing the mean evolutionary tracks of the SSPs (as determined by the route isochrones). The extension must correspond to the 'Z' vector metallicities. Each extension is a data.frame with columns 'SMstar' (fraction of mass in luminous stars) 'SMgas' (fraction of mass in gas) 'SMtot' (fraction of total mass, always 1 in practice) 'SFR' (star formation rate, always 0 in practice) 'SMrem' (fraction of mass in stellar remnants, e.g. black holes etc). Nrow should be equal to the number of age bins (corresponding with 'Age').

Note the "Labels" component of the structure is missing in the FITS version of the spectral library. This is created when processing it with \code{speclib_FITSload}.

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
}
\value{
\code{speclib_download} If \option{stellpop} = 'avail' then the contents of avail.txt will be displayed (using the default method). Otherwise, the target SSP rda/FITS file is downloaded and the destination path (+ filename) is returned.

\code{speclib_FITSload} Returns the target SSP in \code{ProSpect} compatible format. Note this needs to be assigned to a larget object, i.e. \code{newlib = speclib_FITSload('newlib.fits')}.

\code{speclib_check} Returns invisibily TRUE (all checks pass) or FALSE (some checks fail). Also prints the results of various checks to screen.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{ProSpectSED}}
}
\examples{
\dontrun{
# Just check what is there:
speclib_download()

# This is roughly what loading an additional SSP library will look like:
speclib_download(stellpop='XSL_PC_Kr', destpath='~/Downloads/')
load('~/Downloads/XSL_PC_Kr.rda')
temp = ProSpectSED(speclib=XSL_PC_Kr)
plot(temp)

# Check an SSP:
data(BC03lr)
speclib_check(BC03lr)

# Check the FITS version of the library:
BC03lr_file = system.file('BC03lr.fits', package="ProSpectData")
BC03lr_test = speclib_FITSload(BC03lr_file)

#The differences should be very small (parts in 1e-5 - 1e-8):
all.equal(BC03lr, BC03lr_test)
}
}

% \concept{ ~speclib }
% \concept{ ~stellpop }
% \concept{ ~SSP }
% \concept{ ~download }

