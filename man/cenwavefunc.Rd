\name{cenwavefunc}
\alias{cenwavefunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Central Wavelength Calculator
}
\description{
Calculate the central wavelength of a given filter.
}
\usage{
cenwavefunc(filter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
Two column matrix; column one must be flux in Angstroms, column two must be a filter response. These are the type returned by 
\code{\link{getfilt} (see Examples).}
}
}
\details{
The central wavelength (often called the pivot wavelength) is used to compute various filter properties. In ProSpect they have been pre-computed for the provided band passes (see \code{\link{cenwave}}).
}
\value{
Returns the central wavelength in Angstroms.
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cenwave}}
}
\examples{
data('cenwave')
cenwave[8,]
cenwavefunc(getfilt('r_VST'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spectrum }% use one of  RShowDoc("KEYWORDS")
