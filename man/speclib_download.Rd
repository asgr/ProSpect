\name{speclib_download}
\alias{speclib_download}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProSpect Stellar Population Library Downloader
}
\description{
Simple interface to show all the currently available ProSpect SSPs, and download them as needed.
}
\usage{
speclib_download(stellpop = "avail", destpath = "", URL = "https://tinyurl.com/prospect-speclib/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stellpop}{
Character scalar; either 'avail' which will then display all the SSPs currently available at \option{URL}, or one of (currently): 'BC03lr' / 'BC03hr' / 'E'MILES' / 'BPASS' / 'XSL_PC_Kr'. Note only XSL_PC_Kr is not available in the core package, but all are included for completeness (and in case in future some are removed from the core package).
}
  \item{destpath}{
Character scalar; where the target SSP should be downloaded to.
}
  \item{URL}{
Character scalar; location of the ProSpect SSP libraries. This should not need to be altered by the user.
}
}
\details{
This interface has proven necessary in order to provide flexible SSPs without making the core package too big. We have also had issue with the package getting close to the file size limit of GitHub. In future all additional SSP will be provided through this mechanism.
}
\value{
If \option{stellpop} = 'avail' then the contents of avail.txt will be displayed (using the default method). Otherwise, nothing is returned by the side effect is the target SSP rda file is downloaded.
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
speclib_download(stellpop = 'XSL_PC_Kr', destpath = '~/Downloads/')
load('~/Downloads/XSL_PC_Kr.rda')
temp = ProSpectSED(speclib = XSL_PC_Kr)
plot(temp)
}
}

% \concept{ ~speclib }
% \concept{ ~stellpop }
% \concept{ ~SSP }
% \concept{ ~download }

