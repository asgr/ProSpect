\name{magAB2Jansky}
\alias{magAB2Jansky}
\alias{Jansky2magAB}
\alias{magAB2CGS}
\alias{CGS2magAB}
\alias{Jansky2CGS}
\alias{CGS2Jansky}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between AB mag, Jansky and CGS
}
\description{
Simple functions to convert between the common flux and mag systems.
}
\usage{
magAB2Jansky(x)
Jansky2magAB(x)
magAB2CGS(x)
CGS2magAB(x)
Jansky2CGS(x)
CGS2Jansky(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the value to be converted (LHS of the function name).
}
}
\details{
For clarity the conversion used are provided below:

magAB2Jansky=function(x){10^(-0.4*(x-8.9))}
Jansky2magAB=function(x){-2.5*log10(x)+8.9}

magAB2CGS=function(x){10^(-0.4*(x+48.6))}
CGS2magAB=function(x){-2.5*log10(x)-48.6}

Jansky2CGS=function(x){x*1e-23}
CGS2Jansky=function(x){x*1e23}
}
\value{
Numeric vector; the converted value (RHS of the function name).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{photom}}, \code{\link{magABcalc}}
}
\examples{
magAB2Jansky(20)
Jansky2magAB(magAB2Jansky(20))
Jansky2CGS(magAB2Jansky(20))

magAB2CGS(20)
CGS2magAB(magAB2CGS(20))
CGS2Jansky(magAB2CGS(20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ magAB }% use one of  RShowDoc("KEYWORDS")
\keyword{ Jansky }% __ONLY ONE__ keyword per line
\keyword{ CGS }
