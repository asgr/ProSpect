\name{cenwave}
\alias{cenwave}
\alias{pivwave}
\docType{data}
\title{
Filter Central and Pivot Wavelengths
}
\description{
Table of central and pivot wavelengths of all provided filters. Here we use the most common names for the respective filters (e.g. K_VISTA for Ks_VISTA and S500 for 500 micron Spire) and remove duplicates, so be careful when using this table using string logic etc.
}
\usage{
data("cenwave")
data("pivwave")
}
\format{
  A data frame with 39 observations on the following 2 variables.
  \describe{
    \item{\code{filter}}{A character vector}
    \item{\code{cenwave}}{A numeric vector}
  }
}
\examples{
data('cenwave')
print(cenwave)

barplot(height=cenwave$cenwave, names.arg=cenwave$filter, log='y')

\dontrun{
#This is how I make cenwave, assuming we are in the ProSpect main directory.
library(foreach)
filt_files=list.files('data/')
filt_files=filt_files[grep('filt_',filt_files)]
filter=unlist(strsplit(unlist(strsplit(filt_files,'filt_'))[c(F,T)],'.tab'))
cenwave=foreach(i=filt_files,.combine='c')%do%{cenwavefunc(fread(paste0('data/',i)))}

cenwave=data.frame(filter=filter, cenwave=cenwave)
cenwave=cenwave[order(cenwave$cenwave),]
row.names(cenwave)=NULL

save(cenwave, file='data/cenwave.rda')
}
}
\keyword{datasets}
