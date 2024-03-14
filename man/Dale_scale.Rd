\name{Dale_scale}
\alias{Dale_scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale Dale Templates
}
\description{
Correctly rescales Dale templates as needed for given parameters.
}
\usage{
Dale_scale(alpha_SF = 1.5, AGNfrac = 0.5, Dale_in)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha_SF}{
Numeric scalar; desired interpolated alpha slope of the star forming population. Lower values mean hotter dust.
}
  \item{AGNfrac}{
Numeric scalar; AGN fraction linearly mixed in energy over the range 5-20 microns (see Dale 2014).
}
  \item{Dale_in}{
Two column matrix; the wavelength and flux of the dust model to be computed. \option{alpha_SF} and \option{AGNfrac} should still be provided correctly, but passing in \option{Dale_in} will often be faster because there will be less data loading required.
}
}
\details{
Something to be aware of is how Dale 2014 does the mixing. It is actually defined as the energy intergral between 5 and 20 microns. This means the "AGNfrac" as computed is not representative of the total energy in AGN compared to attenuated star formation. If you compare our pure SF normalised to 1 versus pure AGN normalised to 1 over this regime (i.e. Dale_NormTot$Aspec[[1]] versus Dale_NormTot$Aspec[[21]]) these actually need to be further scaled by 0.7545396 (for SF) and 0.2454604 (for AGN) to bring them into full alignment. \code{Dale_scale} does this energy scaling for the user, also see examples.
}
\value{
Returns the two component vector bolometric luminosity scaling for dust (element 1, Dustfrac_bol) and AGN (element 2, AGNfrac_bol) components of the FIR emission.
}
\references{
Dale et al, 2014, ApJ, 784, 11 \cr
Draine et al, 2007, ApJ, 663, 866
}
\author{
Aaron Robotham
}
\seealso{
\code{Dale}, ~~~
}
\examples{
#Check our Dale_scale function for an interpolated AGNfrac and alpha_SF:

alpha_SF = 1.53
AGNfrac = 0.57

Dale_mix = Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, type='NormTot')
Dale_AGN = Dale_interp(alpha_SF=alpha_SF, AGNfrac=1, type='NormTot')
Dale_SFR = Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, type='NormTot')

tempscale = Dale_scale(alpha_SF=alpha_SF, AGNfrac=AGNfrac)

plot(Dale_NormTot$Wave, Dale_mix$Aspec, type='l', log='xy', lwd=4)
lines(Dale_NormTot$Wave, Dale_SFR$Aspec*tempscale[1], col='red')
lines(Dale_NormTot$Wave, Dale_AGN$Aspec*tempscale[2], col='brown')
lines(Dale_NormTot$Wave,
Dale_SFR$Aspec*tempscale[1] + Dale_AGN$Aspec*tempscale[2],
col='green', lty=3, lwd=2)

#The above is how you would reconstruct a generic fit, where you need to scale these
#numbers by the total attenuated flux/luminosity (since above is normalised to 1).
}

