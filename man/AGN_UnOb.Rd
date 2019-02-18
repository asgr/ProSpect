\name{AGN_UnOb}
\alias{AGN_UnOb}
\alias{AGN_UnOb_Sparse}
\docType{data}
\title{
Unattentuated AGN Templates
}
\description{
Unattentuated AGN templates taken from Andrews et al 2018
}
\usage{
data("AGN_UnOb")
data("AGN_UnOb_Sparse")
}
\format{
  AGN_UnOb is a data frame with 6492 observations on the following 2 variables.
  \describe{
    \item{\code{Wave}}{a numeric vector}
    \item{\code{Aspec}}{a numeric vector}
  }
  
  AGN_UnOb_Sparse is a data frame with 650 observations on the following 2 variables.
  \describe{
    \item{\code{Wave}}{a numeric vector}
    \item{\code{Aspec}}{a numeric vector}
  }
  
  In both cases the spectra are normalised so that they integrate to 1 (see examples).
}
\details{
We use the unattentuated AGN template seen in Andrews et al (2018)
}
\source{
Andrews et al, 2018, MNRAS, 474, 898
}
\references{
Andrews et al, 2018, MNRAS, 474, 898
}
\examples{
data(AGN_UnOb)
data(AGN_UnOb_Sparse)
plot(AGN_UnOb, log='xy', type='l')
lines(AGN_UnOb_Sparse, col='red')

sum(c(0,diff(AGN_UnOb$Wave))*AGN_UnOb$Aspec)
sum(c(0,diff(AGN_UnOb_Sparse$Wave))*AGN_UnOb_Sparse$Aspec)
}
\keyword{datasets}
