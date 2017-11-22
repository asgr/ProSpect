\name{CF_general}
\alias{CF_general}
\alias{CF_birth}
\alias{CF_screen}
\alias{CF}
\alias{dust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Charlot and Fall Simple Dust Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CF_general(wave, tau = 0.3, pow = -0.7, pivot=5500)
CF_birth(wave, tau = 1.0,pow = -0.7, pivot=5500)
CF_screen(wave, tau = 0.3, pow = -0.7, pivot=5500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
Observed wavelength in Angstroms.
}
  \item{tau}{
The exponential tau of the dust, where larger means more reddening. Charlot and Fall suggested defaults.
}
  \item{pow}{
Power to further raise the attenuation. Charlot and Fall suggest this should be close to -0.7 based on observations.
}
  \item{pivot}{
The pivot wavelength in Angstroms. Charlot and Fall suggested default.
}
}
\details{
\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (wave, tau = 0.3, pow = -0.7)
{
    return = exp(-tau * (wave/5500)^(pow))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
