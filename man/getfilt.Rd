\name{getfilt}
\alias{getfilt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Requested Filter Curve
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getfilt(filter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filter}{
%%     ~~Describe \code{filter} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function (filter)
{
    out = NA
    if (filter == "FUV") {
        data("filt_FUV_GALEX")
        out = filt_FUV_GALEX
    }
    if (filter == "NUV") {
        data("filt_NUV_GALEX")
        out = filt_NUV_GALEX
    }
    if (filter == "u_SDSS") {
        data("filt_u_SDSS")
        out = filt_u_SDSS
    }
    if (filter == "g_SDSS") {
        data("filt_g_SDSS")
        out = filt_g_SDSS
    }
    if (filter == "r_SDSS") {
        data("filt_r_SDSS")
        out = filt_r_SDSS
    }
    if (filter == "i_SDSS") {
        data("filt_i_SDSS")
        out = filt_i_SDSS
    }
    if (filter == "z_SDSS") {
        data("filt_z_SDSS")
        out = filt_z_SDSS
    }
    if (filter == "u_VST") {
        data("filt_u_VST")
        out = filt_u_VST
    }
    if (filter == "g_VST") {
        data("filt_g_VST")
        out = filt_g_VST
    }
    if (filter == "r_VST") {
        data("filt_r_VST")
        out = filt_r_VST
    }
    if (filter == "i_VST") {
        data("filt_i_VST")
        out = filt_i_VST
    }
    if (filter == "z_VST") {
        data("filt_z_VST")
        out = filt_z_VST
    }
    if (filter == "Z_VISTA") {
        data("filt_Z_VISTA")
        out = filt_Z_VISTA
    }
    if (filter == "Y_VISTA") {
        data("filt_Y_VISTA")
        out = filt_Y_VISTA
    }
    if (filter == "J_VISTA") {
        data("filt_J_VISTA")
        out = filt_J_VISTA
    }
    if (filter == "H_VISTA") {
        data("filt_H_VISTA")
        out = filt_H_VISTA
    }
    if (filter == "K_VISTA") {
        data("filt_K_VISTA")
        out = filt_K_VISTA
    }
    if (filter == "Ks_VISTA") {
        data("filt_Ks_VISTA")
        out = filt_Ks_VISTA
    }
    if (filter == "W1") {
        data("filt_W1_WISE")
        out = filt_W1_WISE
    }
    if (filter == "W2") {
        data("filt_W2_WISE")
        out = filt_W2_WISE
    }
    if (filter == "W3") {
        data("filt_W3_WISE")
        out = filt_W3_WISE
    }
    if (filter == "W4") {
        data("filt_W4_WISE")
        out = filt_W4_WISE
    }
    if (filter == 100 | filter == "100") {
        data("filt_P100_Herschel")
        out = filt_P100_Herschel
    }
    if (filter == 160 | filter == "160") {
        data("filt_P160_Herschel")
        out = filt_P160_Herschel
    }
    if (filter == 250 | filter == "250") {
        data("filt_S250_Herschel")
        out = filt_S250_Herschel
    }
    if (filter == 350 | filter == "350") {
        data("filt_S350_Herschel")
        out = filt_S350_Herschel
    }
    if (filter == 450 | filter == "450" | filter == 500 | filter ==
        "500") {
        data("filt_S450_Herschel")
        out = filt_S450_Herschel
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
