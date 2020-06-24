\name{interp_param}
\alias{interp_param}
\alias{interp_quick}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate weights for a vector of values
}
\description{
Utility functions to interpolate the lineat of log spaced weights for a vector of values. This can be used to create interpolations for x versus y type data, but it is more useful when the property to be weighted is non-trivial, e.g. you have spectra on a grid of metallicities and you want to interpolate between values. \code{interp_quick} is faster, but can only be used for single value lookups, has less flesibility, and returns less information.
}
\usage{
interp_param(x, params, log = FALSE, method = "linear")
interp_quick(x, params, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the value/s to use for the interpolation. For \code{interp_quick} this can only be a scalar.
}
  \item{params}{
Numeric vector; the values to be interpolated.
}
  \item{log}{
Logical scalar; should the interpolation be done in linear space (\option{log}=FALSE, the default), or log space (\option{log}=TRUE).
}
  \item{method}{
One of “constant", “linear", “nearest", “spline", or “cubic"; default is “linear". This is passed to \code{\link{interp1}}.
}
}
\details{
This routine is used to calculate appropriate weights for a number of interpolations within \code{ProSpect}, where outputs are often generated for fixed grids of parameters (e.g. metallicity, AGN fraction and radiation field).
}
\value{
x param_lo param_hi ID_lo weight_lo ID_hi weight_hi flag
A data.frame with the same number of rows as the length of \option{x}, with columns:

\item{x}{The value/s uses for the interpolation (might be different to the input \option{x} if this went beyond the limits of \option{params})}
\item{param_lo}{The nearest value of \option{params} lower than \option{x}}
\item{param_hi}{The nearest value of \option{params} higher than \option{x}}
\item{ID_lo}{The location ID of the nearest value of \option{params} lower than \option{x}}
\item{ID_hi}{The location ID of the nearest value of \option{params} higher than \option{x}}
\item{ID_mode}{The ID with the most weight between \option{ID_lo} and \option{ID_hi}}
\item{wt_lo}{The interpolation weight to be applied to the nearest value of \option{params} lower than \option{x}}
\item{wt_hi}{The interpolation weight to be applied to the nearest value of \option{params} higher than \option{x}}
\item{flag}{Interpolation flag: 0 means the input \option{x} is at an exact value of \option{params}, in which case \option{ID_lo}=\option{ID_hi}, \option{weight_lo}=1 and \option{weight_hi}=0; 1 means \option{x} is less than the minimum of \option{params}, so forced to this value; 2 means \option{x} is between two values in \option{params}, so interpolation behaves in a standard sense; 3 means \option{x} is more than the maximum of \option{params}, so forced to this value}
%% ...
}
\author{
Aaron Robotham
}
\note{
In the output, \option{flag}=2 is the "normal" flag in the sense the interpolation has not gone beyond the limits of \option{params} and is not trivial (as exact value of \option{params}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Dale_interp}}
}
\examples{
interp_param(c(0.1,3.3,5.8,8,11.2),1:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interpolation }% use one of  RShowDoc("KEYWORDS")
