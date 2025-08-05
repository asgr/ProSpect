\name{Inoue14Coefficients}
\alias{Inoue14_LAFcoef}
\alias{Inoue14_DLAcoef}
\docType{data}
\title{
Inoue14 Line Series Coefficients
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Coefficients for the Lyman alpha and damped Lyman series absorption in the Inoue14 model. 
}
\usage{
data("Inoue14_LAFcoef")
data("Inoue14_DLAcoef")
}
\format{

The format of Inoue14_LAFcoef is:\cr
'data.frame':	39 obs. of 5 variables:\cr
   ..$ j   : num [1:39] 2  3  4  5  6 ...\cr
   ..$ lambda: chr [1:47] 1215.670 1025.720  972.537  949.743  937.803 ...\cr
   ..$ AJ1  : chr [1:47] 0.016897600 0.004692290 0.002238980 0.001319010 0.000870656 ...\cr
   ..$ AJ2   : num [1:47] 0.002353790 0.000653625 0.000311884 0.000183735 0.000121280 ...\cr
   ..$ AJ3   : num [1:47] 1.02611e-04 2.84940e-05 1.35962e-05 8.00974e-06 5.28707e-06 ...\cr
   
The format of Inoue14_DLAcoef is:\cr
'data.frame':	39 obs. of 4 variables:\cr
   ..$ j   : num [1:39] 2  3  4  5  6 ...\cr
   ..$ lambda: chr [1:47] 1215.670 1025.720  972.537  949.743  937.803 ...\cr
   ..$ AJ1  : chr [1:47] 0.000161698 0.000154539 0.000149767 0.000146031 0.000142893 ...\cr
   ..$ AJ2   : num [1:47] 5.38995e-05 5.15129e-05 4.99222e-05 4.86769e-05 4.76312e-05 ...\cr
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Inoue (2014) connects the optical depth of the IGM to the distribution of absorbers as a function of redshift and hydrogen column density. The model includes opacity from Lyman continuum and Lyman series line absorption in the Lyman alpha forest (LAF, log(N_H / cm^-2) < 17.2), Lyman limit systems (LLS, 17.2 < log(N_H / cm^-2) <= 20.3) and damped Lyman systems (DLA, log(N_H / cm^-2) > 20.3). 
Eqs. 21-29 in Inoue (2014).
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Inoue14_IGM}}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Data obtained from https://github.com/gbrammer/eazy-py/tree/master/eazy/data
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Inoue, A. K., Shimizu, I., Iwata, I., & Tanaka, M. 2014, Monthly Notices of the Royal Astronomical Society, 442, 1805
}
\examples{
data(Inoue14_LAFcoef)
data(Inoue14_DLAcoef)
}
\keyword{datasets}
