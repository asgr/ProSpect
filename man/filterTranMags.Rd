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
filterTranMags(mag_in, mag_out, return = 'all')
filterTranBands(filt_in, filt_out, zrange = c(0,0.5), Nsamp = 1e3, seed = 666,
return = 'all')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mag_in}{
Numeric matrix, data.frame or vector; columns should be observed magnitudes from bluest (LHS) to reddest (RHS). If \option{mag_in} is a vector then only the offset beta is computed, not colour terms, i.e. \option{mag_out} ~ \option{mag_in} + beta.
}
  \item{mag_out}{
Numeric vector; magnitudes of target filter to map.
}
  \item{return}{
Character scalar; either 'all' (returning all fit combinations); 'bestalpha' (returning just the best in terms of alpha, where smallest alpha term is preferred); 'bestbeta' (returning just the best in terms of beta, where smallest beta term is preferred); 'bestscat' (returning just the best in terms of scatter, where smallest scat term is preferred); 'bestall' (returning the best for all terms, where the smallest sum of 2.alpha^2 + beta^2 + scat^2 is preferred).
}
  \item{filt_in}{
List of functions; filters to be used to calculate filter conversions. In general these should at least be near each other in wavelength. If only one filter function is provided then only the beta offset is computed, i.e. \option{filt_out} ~ \option{filt_in} + beta.
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
\code{filterTranMags} returns the functional mappings for all sensible combinations of \option{mag_in}. This will create a mapping based on a filter mag, colour term and offset. If \option{mag_in} is a vector then only the offset is calculated (it is not possible to compute a colour term), i.e. \option{mag_out} ~ \option{mag_in} + beta. In all cases the scatter around the best solution is also provided (\option{scat}). In general smaller scatter is better, but sometimes it is better to have a smaller \option{alpha} since colour noise can be the dominant error.  \cr

\code{filterTranBands} returns the functional mappings, but it also does the SFH sampling to create the magnitude inputs for \code{filterTranMags}. This means it is the higher level function, and will likely be how most people do the filter transformations.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{getfilt}}, \code{\link{EAZY_filters}}
}
\examples{
\dontrun{
#Using filterTranBands
filters = c('g_SDSS', 'r_SDSS', 'i_SDSS')
filt_in = {}
for(i in filters){filt_in=c(filt_in,list(approxfun(getfilt(i))))}
names(filt_in) = filters
filt_out = approxfun(getfilt('r_VST'))

filterTranBands(filt_in, filt_out)

#Using filterTranMags with Shark SFH samples:

data(ProFiltTrans_Shark)
filterTranMags(ProFiltTrans_Shark$maglist$Age3[,c("g_SDSS", "r_SDSS", "i_SDSS")],
               ProFiltTrans_Shark$maglist$Age3[,"r_VST"])

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ filter }% use one of  RShowDoc("KEYWORDS")
