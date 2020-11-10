\name{LKL10}
\alias{LKL}
\alias{LKL10}
\alias{LKL10_NormHalpha}
\alias{LKL10_NormHbeta}
\alias{LKL10_NormHlines}
\alias{LKL10_NormAll}
\docType{data}
\title{
Select Outputs of Levesque, Kewley & Larson (2010)
}
\description{
Line luminosities as a function of metallicity (\option{Z}) and radiation level (\option{q}) taken from Levesque, Kewley & Larson (2010). These are used by \code{\link{emissionLines}}.
}
\usage{
data("LKL10_NormHalpha")
data("LKL10_NormHbeta")
data("LKL10_NormHlines")
data("LKL10_NormAll")
}
\format{
Each output has the following format internal format:\cr

  The format is:\cr
List of 3\cr
 $ Z       : num [1:5] 0.001 0.004 0.008 0.02 0.04\cr
 $ q       : num [1:7] 1e+07 2e+07 4e+07 8e+07 1e+08 2e+08 4e+08\cr
 $ Mappings:List of 5\cr
  ..$ Z0.001:'data.frame':	47 obs. of  10 variables:\cr
  .. ..$ wave   : num [1:47] 1216 2325 2798 3726 3729 ...\cr
  .. ..$ element: chr [1:47] "H" "C" "Mg" "O" ...\cr
  .. ..$ state  : chr [1:47] "I" "II" "II" "II" ...\cr
  .. ..$ q1e7   : num [1:47] 11.03806 0.00869 0.01477 0.08713 0.11765 ...\cr
  .. ..$ q2e7   : num [1:47] 10.91982 0.00941 0.01306 0.09997 0.13481 ...\cr
  .. ..$ q4e7   : num [1:47] 10.76628 0.00928 0.01064 0.10118 0.13629 ...\cr
  .. ..$ q8e7   : num [1:47] 10.28846 0.00771 0.00751 0.08101 0.10896 ...\cr
  .. ..$ q1e8   : num [1:47] 10.24399 0.00736 0.00684 0.07546 0.10144 ...\cr
  .. ..$ q2e8   : num [1:47] 10.11088 0.00614 0.00504 0.05696 0.07626 ...\cr
  .. ..$ q4e8   : num [1:47] 9.93377 0.00485 0.00361 0.04008 0.05329 ...\cr
  ..$ Z0.004:'data.frame':	47 obs. of  10 variables:\cr
  .. ..$ wave   : num [1:47] 1216 2325 2798 3726 3729 ...\cr
  .. ..$ element: chr [1:47] "H" "C" "Mg" "O" ...\cr
  .. ..$ state  : chr [1:47] "I" "II" "II" "II" ...\cr
  .. ..$ q1e7   : num [1:47] 9.3835 0.0404 0.0437 0.2744 0.3706 ...\cr
  .. ..$ q2e7   : num [1:47] 9.0943 0.0446 0.0361 0.3066 0.4137 ...\cr
  .. ..$ q4e7   : num [1:47] 8.8008 0.0448 0.0262 0.3188 0.4302 ...\cr
  .. ..$ q8e7   : num [1:47] 8.5762 0.0417 0.0177 0.3082 0.416 ...\cr
  .. ..$ q1e8   : num [1:47] 8.496 0.0405 0.0155 0.3025 0.4081 ...\cr
  .. ..$ q2e8   : num [1:47] 8.4438 0.0376 0.0107 0.2854 0.3843 ...\cr
  .. ..$ q4e8   : num [1:47] 8.50808 0.03764 0.00839 0.28749 0.38398 ...\cr
  ..$ Z0.008:'data.frame':	47 obs. of  10 variables:\cr
  .. ..$ wave   : num [1:47] 1216 2325 2798 3726 3729 ...\cr
  .. ..$ element: chr [1:47] "H" "C" "Mg" "O" ...\cr
  .. ..$ state  : chr [1:47] "I" "II" "II" "II" ...\cr
  .. ..$ q1e7   : num [1:47] 8.8061 0.0683 0.0633 0.3985 0.5385 ...\cr
  .. ..$ q2e7   : num [1:47] 8.6462 0.067 0.0516 0.4089 0.5523 ...\cr
  .. ..$ q4e7   : num [1:47] 8.4942 0.0586 0.0375 0.3698 0.4995 ...\cr
  .. ..$ q8e7   : num [1:47] 8.4137 0.0482 0.0259 0.3044 0.4102 ...\cr
  .. ..$ q1e8   : num [1:47] 8.3986 0.0449 0.0228 0.2819 0.3794 ...\cr
  .. ..$ q2e8   : num [1:47] 8.3735 0.036 0.0159 0.2194 0.2927 ...\cr
  .. ..$ q4e8   : num [1:47] 8.4207 0.0299 0.0122 0.1784 0.2324 ...\cr
  ..$ Z0.02 :'data.frame':	47 obs. of  10 variables:\cr
  .. ..$ wave   : num [1:47] 1216 2325 2798 3726 3729 ...\cr
  .. ..$ element: chr [1:47] "H" "C" "Mg" "O" ...\cr
  .. ..$ state  : chr [1:47] "I" "II" "II" "II" ...\cr
  .. ..$ q1e7   : num [1:47] 8.101 0.036 0.026 0.231 0.311 ...\cr
  .. ..$ q2e7   : num [1:47] 8.1007 0.0331 0.0202 0.2249 0.3024 ...\cr
  .. ..$ q4e7   : num [1:47] 8.0986 0.0293 0.0146 0.2114 0.284 ...\cr
  .. ..$ q8e7   : num [1:47] 8.1033 0.0258 0.0104 0.1936 0.2595 ...\cr
  .. ..$ q1e8   : num [1:47] 8.10398 0.02486 0.00935 0.18753 0.251 ...\cr
  .. ..$ q2e8   : num [1:47] 8.146 0.0232 0.0073 0.1742 0.2293 ...\cr
  .. ..$ q4e8   : num [1:47] 8.24788 0.02345 0.00662 0.1727 0.21722 ...\cr
  ..$ Z0.04 :'data.frame':	47 obs. of  10 variables:\cr
  .. ..$ wave   : num [1:47] 1216 2325 2798 3726 3729 ...\cr
  .. ..$ element: chr [1:47] "H" "C" "Mg" "O" ...\cr
  .. ..$ state  : chr [1:47] "I" "II" "II" "II" ...\cr
  .. ..$ q1e7   : num [1:47] 7.9933 0.0173 0.0131 0.1306 0.1751 ...\cr
  .. ..$ q2e7   : num [1:47] 7.98401 0.01236 0.00887 0.1046 0.14035 ...\cr
  .. ..$ q4e7   : num [1:47] 7.97742 0.00846 0.00564 0.07615 0.10213 ...\cr
  .. ..$ q8e7   : num [1:47] 7.97422 0.00615 0.00373 0.05516 0.0736 ...\cr
  .. ..$ q1e8   : num [1:47] 7.9835 0.00582 0.00343 0.05149 0.06832 ...\cr
  .. ..$ q2e8   : num [1:47] 8.05199 0.00562 0.00309 0.04768 0.06094 ...\cr
  .. ..$ q4e8   : num [1:47] 8.21028 0.00651 0.00348 0.05369 0.06299 ...\cr
}
\details{
Select outputs of Levesque, Kewley & Larson (2010). Here we selected an electron density of 100. Sub lists are sorted by metallicity (as per \option{Z}) and within each list are columns specifying the wavelength (Ang), element, state (neutral line = I, single ionised = II, double ionised = III etc), and radiation field as per \option{q}.

There are 4 variants each scales differently. \code{LKL10_NormHalpha} is scale so that the Halpha line is 1, and all other lines are scaled in proportion to this. \code{LKL10_NormHbeta} is scale so that the Hbeta line is 1, and all other lines are scaled in proportion to this. \code{LKL10_NormHlines} is scale so that the sum of all hydrogen lines is 1, and all other lines are scaled in proportion to this. \code{LKL10_NormAll} is scale so that the sum of lines is 1.
}
\source{
Levesque, Kewley & Larson, 2010, AJ, 139, 712
}
\references{
Levesque, Kewley & Larson, 2010, AJ, 139, 712
}
\examples{
data(LKL10_NormHalpha)
data(LKL10_NormHbeta)
data(LKL10_NormHlines)
data(LKL10_NormAll)
}
\keyword{datasets}
