\name{SFHp4}
\alias{SFHp4}
\alias{SMstarp4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Four Phase Star Formation History
}
\description{
Computes output galaxy spectra using four distinct phases of constant star formation: recent burst, young, old and ancient. Each can have a different metallicity associated with it (default solar).
}
\usage{
SFHp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop = "BC03lr", tau_birth = 1, tau_screen = 0.3,
filters = "all", Z = c(5, 5, 5, 5), z = 0.1, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = "planck",
outtype = "mag", cossplit = c(9e+09, 1.3e+10), dosplit = FALSE)

SMstarp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), Z = c(5, 5, 5, 5), z = 0, H0 = 100, OmegaM = 0.3,
OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = 'planck',
cossplit = c(9e9, 1.3e10), dosplit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burstmass}{
Numeric scalar; the total stellar mass formed in a recent burst.
}
  \item{youngmass}{
Numeric scalar; the total stellar mass formed in young stars.
}
  \item{oldmass}{
Numeric scalar; the total stellar mass formed in old stars.
}
  \item{ancientmass}{
Numeric scalar; the total stellar mass formed in ancient stars.
}
  \item{burstage}{
Numeric vector; the age limits in years for the burst phase.
}
  \item{youngage}{
Numeric vector; the age limits in years for the young phase.
}
  \item{oldage}{
Numeric vector; the age limits in years for the old phase.
}
  \item{ancientage}{
Numeric vector; the age limits in years for the ancient phase.
}
  \item{tau_birth}{
Numeric scalar; dust tau for birth clouds. Associated with sub 100 Myr star formation.
}
  \item{tau_screen}{
Numeric scalar; dust tau for screen.
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters.
}
  \item{Z}{
Integer vector of length 4; specifies which metallicity to use for each star formation phase.
}
  \item{z}{
Numeric scalar; cosmological redshift, where z must be > -1.
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=100 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.3).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1-OmegaM-OmegaR = 0.7).
}
  \item{OmegaR}{
Numeric scalar; Omega Radiation today (default is 0, but OmegaM/3400 is typical).
}
  \item{w0}{
The value of dark energy equation of state at z=0. See \code{\link{cosgrow}} for more details.
}
  \item{wprime}{
The evolution term that governs how the dark energy equation of state evolves with redshift. See \code{\link{cosgrow}} for more details.
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck=Planck15 and WMAP=WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
}
  \item{outtype}{
Character scalar; the type of output to generate. Options are mag, cgs and janksky. See \code{\link{magABcalc}} for details.
}
  \item{cossplit}{
Numeric vector; the cosmological age at which the ancient star formation event begins and ends (only used if \option{dosplit=TRUE}).
}
  \item{dosplit}{
Boolean scalar; should the \option{cossplit} ancient star formation limits be used?
}
}
\details{
Optionally the definition of old and ancient can be associated with a fixed time the history of the Universe (e.g. the peak of starformation). This means the exact lookback period adapts as a function of the redshift provided. This often makes more sense then simply fixing the age of the ancient population since it precents the stellar populations becoming older than the Universe.
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
Aaron Robotham
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
function (burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, 
    ancientmass = 1e+10, burstage = c(0, 1e+08), youngage = c(1e+08, 
        1e+09), oldage = c(1e+09, 9e+09), ancientage = c(9e+09, 
        1.3e+10), tau_birth = 1, tau_screen = 0.3, filters = "all", 
    Z = c(5, 5, 5, 5), z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - 
        OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref = "planck", 
    outtype = "mag", cossplit = c(9e+09, 1.3e+10), dosplit = FALSE) 
{
    if (length(Z) == 1) {
        Z = rep(Z, 4)
    }
    if (dosplit) {
        TravelTime = cosdistTravelTime(z = z, ref = "planck") * 
            1e+09
        Tsplit = cossplit[1] - TravelTime
        Tstart = cossplit[2] - TravelTime
        oldage[2] = Tsplit
        ancientage[1] = Tsplit
        ancientage[2] = Tstart
    }
    speclib_burst = BC03lr$Zspec[[Z[1]]]
    speclib_young = BC03lr$Zspec[[Z[2]]]
    speclib_old = BC03lr$Zspec[[Z[3]]]
    speclib_ancient = BC03lr$Zspec[[Z[4]]]
    if (tau_birth != 0) {
        speclib_burst[1:70, ] = t(t(speclib_burst[1:70, ]) * 
            CF_birth(BC03lr$Wave, tau = tau_birth))
        speclib_young[1:70, ] = t(t(speclib_young[1:70, ]) * 
            CF_birth(BC03lr$Wave, tau = tau_birth))
        speclib_old[1:70, ] = t(t(speclib_old[1:70, ]) * CF_birth(BC03lr$Wave, 
            tau = tau_birth))
        speclib_ancient[1:70, ] = t(t(speclib_ancient[1:70, ]) * 
            CF_birth(BC03lr$Wave, tau = tau_birth))
    }
    if (tau_screen != 0) {
        speclib_burst = t(t(speclib_burst) * CF_screen(BC03lr$Wave, 
            tau = tau_screen))
        speclib_young = t(t(speclib_young) * CF_screen(BC03lr$Wave, 
            tau = tau_screen))
        speclib_old = t(t(speclib_old) * CF_screen(BC03lr$Wave, 
            tau = tau_screen))
        speclib_ancient = t(t(speclib_ancient) * CF_screen(BC03lr$Wave, 
            tau = tau_screen))
    }
    burstageloc = c(which.min(abs(BC03lr$Age - burstage[1])), 
        which.min(abs(BC03lr$Age - burstage[2])))
    youngageloc = c(which.min(abs(BC03lr$Age - youngage[1])), 
        which.min(abs(BC03lr$Age - youngage[2])))
    oldageloc = c(which.min(abs(BC03lr$Age - oldage[1])), which.min(abs(BC03lr$Age - 
        oldage[2])))
    ancientageloc = c(which.min(abs(BC03lr$Age - ancientage[1])), 
        which.min(abs(BC03lr$Age - ancientage[2])))
    burstlum = colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2], 
        ]) * BC03lr$AgeWeights[burstageloc[1]:burstageloc[2]]) * 
        burstmass/sum(BC03lr$AgeWeights[burstageloc[1]:burstageloc[2]])
    younglum = colSums(rbind(speclib_young[youngageloc[1]:youngageloc[2], 
        ]) * BC03lr$AgeWeights[youngageloc[1]:youngageloc[2]]) * 
        youngmass/sum(BC03lr$AgeWeights[youngageloc[1]:youngageloc[2]])
    oldlum = colSums(rbind(speclib_old[oldageloc[1]:oldageloc[2], 
        ]) * BC03lr$AgeWeights[oldageloc[1]:oldageloc[2]]) * 
        oldmass/sum(BC03lr$AgeWeights[oldageloc[1]:oldageloc[2]])
    ancientlum = colSums(rbind(speclib_ancient[ancientageloc[1]:ancientageloc[2], 
        ]) * BC03lr$AgeWeights[ancientageloc[1]:ancientageloc[2]]) * 
        ancientmass/sum(BC03lr$AgeWeights[ancientageloc[1]:ancientageloc[2]])
    lum = burstlum + younglum + oldlum + ancientlum
    lumtot = sum(c(0, diff(BC03lr$Wave)) * lum)
    flux = Lum2Flux(wave = BC03lr$Wave, lum = lum, z = z, H0 = H0, 
        OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, 
        wprime = wprime, ref = ref)
    if (filters[1] == "all") {
        filters = cenwave$filter
    }
    out = {
    }
    if (outtype == "mag") {
        for (i in filters) {
            out = c(out, magABcalc(flux, filter = i))
        }
    }
    if (outtype == "cgs" | outtype == "CGS") {
        for (i in filters) {
            out = c(out, CGScalc(flux, filter = i))
        }
    }
    if (outtype == "jansky" | outtype == "Jansky") {
        for (i in filters) {
            out = c(out, Janskycalc(flux, filter = i))
        }
    }
    out = cbind(cenwave[match(filters, cenwave$filter), ], out = out)
    masstot = burstmass + youngmass + oldmass + ancientmass
    ages = rbind(c(burstage, diff(burstage), mean(burstage)), 
        c(youngage, diff(youngage), mean(youngage)), c(oldage, 
            diff(oldage), mean(oldage)), c(ancientage, diff(ancientage), 
            mean(ancientage)))
    colnames(ages) = c("lo", "hi", "duration", "mean")
    masses = cbind(c(burstmass, youngmass, oldmass, ancientmass), 
        c(burstmass + youngmass + oldmass + ancientmass, youngmass + 
            oldmass + ancientmass, oldmass + ancientmass, ancientmass))
    colnames(masses) = c("Forming", "Formed")
    SFR = masses[, "Forming"]/ages[, "duration"]
    sSFR = SFR/masses[, "Formed"]
    return = list(flux = flux, lum = lum, out = out, masstot = masstot, 
        lumtot = lumtot, M2L = masstot/lumtot, call = match.call(), 
        ages = ages, masses = masses, SFR = SFR, sSFR = sSFR)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
