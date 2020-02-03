\name{filterTran}
\alias{filterTranMags}
\alias{filterTranBands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter Transformation Functions
}
\description{
High and low level interface to achieving colour and simple offset filter conversions. These functions can work on simulated or observed data.
}
\usage{
filterTranMags(mag_in, mag_out, return='all')
filterTranBands(filt_in, filt_out, zrange=c(0,0.5), Nsamp=1e3, seed=666)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mag_in}{
Numeric matrix, data.frame or vector; columns should be observed magnitudes from bluest (LHS) to reddest (RHS).
}
  \item{mag_out}{
Numeric vector; magnitudes of target filter to map.
}
  \item{return}{
Character scalar; either 'all' (returning all fit combinations), or 'best' (returning just the best in terms of scatter).
}
  \item{filt_in}{
List of functions; filters to be used to calculate filter conversions. In general these should at least be near each other in wavelength.  
}
  \item{filt_out}{
Function; filter to be targetted for filter conversions.
}
  \item{zrange}{
Numeric vector; length=2, specifying the lower and upper limits of the redshift range to sample.
}
  \item{Nsamp}{
Integer scalar; number of sample to make.  
}
  \item{seed}{
Integer scalar; random seed to use for sampling.  
}
}
\details{
These functions are used to find filter conversions.
}
\value{
\code{filterTranMags} returns the functional mappings for all sensible combinations of \option{mag_in}. This will create a mapping based on a filter mag, colour term and offset. If \option{mag_in} is a vector then only the offset is calculated (it is not possible to compute a colour term). In all cases the scatter around the best solution is also provided (\option{scat}). In general smaller scatter is better, but sometimes it is better to have a smaller \option{alpha} since colour noise can be the dominant error.  \cr

\code{filterTranBands} returns the funcational mappings, but it also does the SFH sampling to create the magnitude inputs for \code{filterTranMags}. This means it is the higher level function, and will likely be how most people do the filter transformations.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{getfilt}}, \code{\link{EAZY_filters}}
}
\examples{
\dontrun{
filters = c('g_SDSS', 'r_SDSS', 'i_SDSS')
filt_in = {}
for(i in filters){filt_in=c(filt_in,list(approxfun(getfilt(i))))}
names(filt_in) = filters
filt_out = approxfun(getfilt('r_VST'))

filterTranBands(filt_in, filt_out)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line