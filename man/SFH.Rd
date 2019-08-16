\name{Star Formation}
\alias{SFH}
\alias{SFHburst}
\alias{SFHp4}
\alias{SFHp5}
\alias{SFHfunc}
\alias{SM}
\alias{SMstarp4}
\alias{SMstarp5}
\alias{SMstarfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate SED from SFH
}
\description{
As of v0.8.0 \code{SFHp4} and \code{SFHp5} are deprecated and will be removed entirely in future versions (probably v1.0.0+). Please use \code{\link{massfunc_b5}} in combination with \code{SFHfunc} (this is now default) to create stepped SFHs. The mass formed in bins is now specified as a SFR, not as a total mass formed.

The star formation history functions discussed here compute output galaxy spectra  and stellar masses using variable star formation histories. This can be either a pure burst (\code{SFHburst}), distinct phases (\code{SFHp4} / \code{SFHp5}), or a functional form (\code{SFHfunc}). Uses either BC03 or EMILES libraries. For \code{SFHp4} and \code{SFHp5}, each star formation phase can have a different metallicity associated with it (default solar), whereas \code{SFHfunc} can take a functional form for the metallicity.
}
\usage{
SFHfunc(massfunc = massfunc_b5, forcemass = FALSE, agescale = 1, stellpop = "BC03lr",
speclib = NULL, tau_birth = 1, tau_screen = 0.3, pow_birth = -0.7, pow_screen = -0.7,
filters = "all", Z = 5, emission = FALSE, veldisp = 50, emission_scale = "FUV",
escape_frac = 1 - emission, Ly_limit = 911.8, LKL10 = NULL, z = 0.1, H0 = 67.8,
OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype = "mag", sparse = 5, intSFR = FALSE,
unimax = 1.38e+10, agemax = NULL, ...) 

SMstarfunc(massfunc = massfunc_b5, forcemass = FALSE, agescale = 1,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), midage = c(1e+09, 5e+09),
oldage = c(5e+09, 9e+09), ancientage = c(9e+09, 1.3e+10), stellpop = "BC03lr",
speclib = NULL, Z = 5, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref,
unimax = 1.38e+10, agemax = NULL, ...)

SFHburst(burstmass = 1e8, burstage = 0, stellpop = 'BC03lr', speclib = NULL,
tau_birth = 1.0, tau_screen = 0.3, pow_birth = -0.7, pow_screen = -0.7,  filters = 'all',
Z = 0.02, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype = 'mag',
sparse = 5, unimax=  13.8e9, agemax = NULL, ...)

SFHp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop = "BC03lr", speclib = NULL, tau_birth = 1,
tau_screen = 0.3, pow_birth = -0.7, pow_screen = -0.7, filters = "all", Z = c(5, 5, 5, 5),
z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype = "mag",
cossplit = c(9e+09, 1.3e+10), dosplit = FALSE, sparse = 5, unimax = 1.38e+10,
agemax = NULL, ...)

SFHp5(burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10, oldmass = 1e+10,
ancientmass = 1e+10, burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09),
midage = c(1e+09, 5e+09), oldage = c(5e+09, 9e+09), ancientage = c(9e+09, 1.3e+10),
stellpop = "BC03lr", speclib = NULL, tau_birth = 1, tau_screen = 0.3, pow_birth = -0.7,
pow_screen = -0.7, filters = "all", Z = c(5, 5, 5, 5, 5), z = 0.1, H0 = 67.8,
OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype = "mag", cossplit = c(9e+09, 1.3e+10),
dosplit = FALSE, sparse = 5, unimax = 1.38e+10, agemax = NULL, ...)

SMstarp4(burstmass = 1e+08, youngmass = 1e+09, oldmass = 1e+10, ancientmass = 1e+10,
burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09), oldage = c(1e+09, 9e+09),
ancientage = c(9e+09, 1.3e+10), stellpop='BC03lr', speclib = NULL, Z = c(5, 5, 5, 5),
z = 0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit = c(9e+09, 1.3e+10),
dosplit = FALSE, unimax = 1.38e+10, agemax = NULL, ...)

SMstarp5(burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10, oldmass = 1e+10,
ancientmass = 1e+10, burstage = c(0, 1e+08), youngage = c(1e+08, 1e+09),
midage = c(1e+09, 5e+09), oldage = c(5e+09, 9e+09), ancientage = c(9e+09, 1.3e+10),
stellpop = "BC03lr", speclib = NULL, Z = c(5, 5, 5, 5, 5), z = 0, H0 = 67.8,
OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit = c(9e+09, 1.3e+10), dosplit = FALSE,
unimax = 1.38e+10, agemax = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burstmass}{
Numeric scalar; the total stellar mass formed in a recent burst in units Msol.
}
  \item{youngmass}{
Numeric scalar; the total stellar mass formed in young stars in units Msol.
}
  \item{midmass}{
Numeric scalar; the total stellar mass formed in mid stars in units Msol.
}
  \item{oldmass}{
Numeric scalar; the total stellar mass formed in old stars in units Msol.
}
  \item{ancientmass}{
Numeric scalar; the total stellar mass formed in ancient stars in units Msol.
}
  \item{massfunc}{
Function; a function which will take \option{age} (in years) as an input and return SFR (in Msun/yr). Additional \option{massfunc} function arguments can be passed down via the \dots (useful for fitting), so be careful to give additional arguments unique names compared to \code{SFHfunc}. Various sensible functions are available in \code{ProSpect}, see \code{\link{massfunc}}.
}
  \item{burstage}{
Numeric vector; the age limits in years for the burst phase (\code{SFHp4} / \code{SFHp5}). Numeric scalar; the age of the burst in years (\code{SFHburst}).
}
  \item{youngage}{
Numeric vector; the age limits in years for the young phase.
}
  \item{midage}{
Numeric vector; the age limits in years for the mid phase.
}
  \item{oldage}{
Numeric vector; the age limits in years for the old phase.
}
  \item{ancientage}{
Numeric vector; the age limits in years for the ancient phase.
}
  \item{forcemass}{
Logical or numeric scalar; if FALSE then the stellar mass is directly computed from the \option{massfunc} function, i.e. the output of this must be a true star formation rate in Msun per year. If \option{forcemass} is a number then this specifies the required integrated stellar mass that must be formed over all time, and the output of \option{massfunc} will linearly scaled as required to achieve this.
}
  \item{agescale}{
Numeric scalar; scaling to apply to the ages (which will passed in years) to \option{massfunc}. E.g. if left at \option{agescale}=1 then \option{massfunc} is expecting the age in years, but if set to 1e-6/1e-9 then \option{massfunc} is expecting the age in Myrs/Gyrs.
}
  \item{stellpop}{
Character scalar; which stellar population library to use, one of \code{\link{BC03lr}} (default), \code{\link{BC03hr}}, \code{\link{EMILES}}.
}
  \item{speclib}{
Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Doing this speeds up the compute time, since there is no need to lazy load from the package.
}
  \item{tau_birth}{
Numeric scalar; dust tau of birth clouds. Associated with sub 10 Myr star formation. See \code{\link{CF_birth}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF_screen}}.
}
  \item{pow_birth}{
Numeric scalar; power to further raise the dust attenuation of birth clouds. See \code{\link{CF_birth}}.
}
  \item{pow_screen}{
Numeric scalar; power to further raise the dust attenuation of the dust screen. See \code{\link{CF_screen}}.
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters. Can also be a list of matrices or functions to pass into \code{\link{bandpass}}. This creates increasingly faster code, but requires more user effort.
}
  \item{Z}{
Numeric scalar (\code{SFHburst}), or integer vector of length 4/5 (\code{SFHp4} / \code{SFHp5}), or a function (\code{SFHfunc}); specifies which metallicity to use for each star formation phase. The metallicity is the vector location \option{Z} in BC03lr$Z / BC03hr$Z / EMILES$Z., i.e. the default of 5 for BC03lr is close to solar (BC03lr$Z[5]=0.02). For \code{SFHfunc} the functional form takes age as an input (in units of \option{agescale}) and outputs the target metallicity in terms of Z (so 0.02 for solar). \code{SFHp4}) and \code{SFHp5} also take functional inputs, in which case the nearest matching discrete metallicity is used for each phase. Additional \option{Z} function arguments can be passed down via the \dots (useful for fitting), so be careful to give additional arguments unique names compared to \code{SFHp4} / \code{SFHp5} /\code{SFHfunc}. For convenience, if a non-integer value is provided then this is interpreted as the desired metallicity value (so Z=0.02 for solar etc), and converted into the appropriate function internally, e.g. \option{Z}=5 and \option{Z}=0.02 will give the same outputs for \code{\link{BC03}}.
}
  \item{emission}{
Logical; should emission features be added to the spectrum? This will be done using \code{\link{emissionLines}}. 
}
  \item{veldisp}{
Numeric scalar; velocity dispersion (km/s) of emission features. See \code{\link{emissionLines}}. 
}
  \item{emission_scale}{
Character; one of 'FUV' or 'SFR'. 'FUV' scales the emission features via an energy balance, absorbing all unattenuated flux below \option{Ly_limit} by a factor (1-\option{escape_frac}) and re-emitting as lines. 'SFR' scales Halpha from Kennicutt (1998), but still absorbs flux below \option{Ly_limit} by a factor (1-\option{escape_frac}). If providing \option{LKL10} to speed code up then it must be \option{LKL10}=\code{\link{LKL10_NormAll}} when \option{emission_scale}='FUV', and \option{LKL10}=\code{\link{LKL10_NormHalpha}} when \option{emission_scale}='SFR'. Otherwise these are loaded internally as appropriate.
}
  \item{escape_frac}{
Numeric scalar; the fraction of photons allows to escape when the wavelength is shorter than \option{Ly_limit}. If \option{emission}=TRUE and \option{emission_scale}='FUV' then the absorbed flux (1-\option{escape_frac}) is re-emitted as lines. By default it is \option{escape_frac} = 1 - \option{emission}, i.e. if we do not want emission then we allow all UV photons out of the star forming region (to possibly then be attenuated by dust), but if we want emission then all of these photons are used to excite our emission lines.
}
  \item{Ly_limit}{
Numeric scalar; the Lyman limit to use (Angstroms). Photons with wavelengths shorted than \option{Ly_limit} are absorbed by a factor 1-\option{escape_frac}. If \option{emission}=TRUE and \option{emission_scale}='FUV' then the absorbed flux (1-\option{escape_frac}) is re-emitted as lines. The \option{Ly_limit} should generally be 911.8 (the actual ionisation limit n=Inf of 13.6 eV), or sometimes 1215.67 (Lyman alpha n=2).
}
  \item{LKL10}{
List; LKL template to use. If \option{LKL10}=NULL then the appropriate template is loaded internally (this is the safest option). To speed things up the user can pass in the appropriate \code{\link{LKL10}} template, but care must be taken to ensure the scaling used makes use of the correct dataset. See \code{\link{emissionLines}}.  
}
  \item{z}{
Numeric scalar; cosmological redshift of the source. If the provided \option{z}=0 then the returned fluxes will be absolute, i.e. absolute magnitudes if \option{outtype}='mag'. If \option{z}<0 then only luminosities are computed (flux parts are skipped).
}
  \item{H0}{
Numeric scalar; Hubble constant as defined at z=0 (default is H0=67.8 (km/s)/Mpc).
}
  \item{OmegaM}{
Numeric scalar; Omega Matter today (default is 0.308).
}
  \item{OmegaL}{
Numeric scalar; Omega Lambda today (default is for a flat Universe with OmegaL = 1 - OmegaM = 0.692).
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / Planck18 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck = Planck18 and WMAP = WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
}
  \item{unimax}{
Numeric scalar; maximum allowed age of any stellar population relative to z=0 (i.e. today). Any star formation that is older than this will be set to 0. Overridden by \option{agemax} if that is provided.
}
  \item{agemax}{
Numeric scalar; maximum allowed age of any stellar population relative to the redshift of the observed galaxy. Any star formation that is older than this will be set to 0. Overridden by \option{unimax} if provided. It is computational a bit faster to provide this since we do not need to compute a lookup time from the redshift.
}
  \item{outtype}{
Character scalar; specifies the type of photometry computed, must be one of 'mag' / 'magAB' (for AB magnitudes, the default), 'Jansky' / 'Jy' (for Jansky fluxes) or 'CGS' (for CGS fluxes). See \code{\link{magABcalc}} for details.
}
  \item{cossplit}{
Numeric vector; the cosmological age at which the ancient star formation event begins and ends (only used if \option{dosplit=TRUE}).
}
  \item{dosplit}{
Boolean scalar; should the \option{cossplit} ancient star formation limits be used?
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make. Higher values mean coarser spectra is processed, which means less accurate photometry, but faster processing. Generally done when computing higher resolution libraries for AB mags.
}
  \item{intSFR}{
Boolean scalar; should the \option{massfunc} be intergrated between ages? This might be necessary if the SFH is quite bursty, but it is more expensive to compute and not required if the SFH is quite smooth.  
}
  \item{\dots}{
Further arguments to be passed to \code{massfunc} or \code{Z} functions, as appropriate (this is achieved by formals matching, so be careful to used unique names).
}
}
\details{
Optionally the definition of old and ancient can be associated with a fixed time the history of the Universe (e.g. the peak of starformation). This means the exact lookback period adapts as a function of the redshift provided. This often makes more sense then simply fixing the age of the ancient population since it precents the stellar populations becoming older than the Universe.

We use the simple dust models for birth clouds and dust screens as given by Charlot and Fall (2000):

\deqn{A=\exp(-\tau*(\lambda/\lambda_{piv})^n)}{A=exp(-tau*(wave/pivot)^pow)}

The defaults should be reasonable in the regime they are used over, where birth cloud dust should only be applied to stellar populations younger than 10 Myrs. Note that BC03 has 70 stellar population spectra which are younger than this, whilst EMILES has none. This means that EMILES cannot realistically capture the birth dust attenuation phase, so for highly star bursting populations BC03 should probably be preferred.
}
\value{
\code{SFHburst} / \code{SFHp4} / \code{SFHp5} / \code{SFHfunc} returns a list containing:

\item{flux}{Numeric Matrix; first column is observed wavelength (Angstroms) and second column is observed flux (erg/s/cm^2/Angstrom).}
\item{out}{Data frame; first column is the filter name, second columns is the central wavelength, third column is the requested output (e.g. AB mag).}
\item{wave_lum}{Numeric vector; luminosity wavelength in Angstroms.}
\item{lum_atten}{Numeric vector; attenuated stellar luminosity in the standard BC03 units of Lsol / Angstrom.}
\item{lum_unatten}{Numeric vector; unattenuated stellar luminosity in the standard BC03 units of Lsol / Angstrom.}
\item{lumtot_unatten}{Numeric scalar, total unattenuated stellar luminosity in Lsol.}
\item{lumtot_atten}{Numeric scalar, total attenuated stellar luminosity in Lsol.}
\item{lumtot_birth}{Numeric scalar, total attenuated luminosity in birth cloud in Lsol.}
\item{lumtot_screen}{Numeric scalar, total attenuated luminosity in screen cloud in Lsol.}
\item{masstot}{Numeric scalar; total mass formed in Msol.}
\item{massvec}{Numeric vector; the total mass formed in each age window  (\code{SFHfunc} only).}
\item{masses}{Numeric Matrix; mass info where rows are the 4/5 phases of star formation, columns are the total formed in that phase [1] and the total formed over all time [2] (\code{SFHp4} / \code{SFHp5} only).}
\item{M2L}{Numeric scalar; total mass to light ratio.}
\item{ages}{Numeric Matrix; age info where rows are the four phases of star formation, columns are the start [1], end [2], duration [3] and mean time [4] of the phases in years (\code{SFHp4} / \code{SFHp5} only). For \code{SFHburst} it is just the input \option{burstage}.}
\item{SFR}{Numeric vector; the star formation rate averaged out over each star formation phase (\code{SFHp4} / \code{SFHp5} only).}
\item{sSFR}{Numeric vector; the specific star formation rate averaged out over each star formation phase (\code{SFHp4} / \code{SFHp5} only).}
\item{SFRburst}{Numeric scalar; the star formation rate over the last 10^8 years.}
\item{emission_input}{List; the inputs sent to \code{\link{emissionLines}} when \option{emission}=TRUE.}

\code{SMstarp4} / \code{SMstarp5} / \code{SMstarfunc} produces both the stellar mass formed (4/5/5 elements) and the stellar mass remaining (a further 4/5/5 elements) for the 5 different time phases provided in the arguments. Finally it produces the total stellar mass formed, and the total stellar mass remaining (not formed- this is the last element listed).
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{ProSpectSED}}, \code{\link{photom}}, \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}, \code{\link{massfunc}}, \code{\link{emissionLines}}
}
\examples{
#A pretty full example showing some low level things:

SFHdemo_dust=SFHp4()
plot(SFHdemo_dust$out[,2:3], log='x', ylim=c(26,17), xlab=BC03lr$Labels$Wavelab,
ylab='Mag', col=rev(rainbow(26, end=2/3)), pch=16)

SMdemo=SMstarp4()
SMdemo

SFHdemo_nodust=SFHp4(tau_birth=0, tau_screen=0)
wave_star=SFHdemo_nodust$flux[,1]

total_atten=sum(c(0,diff(wave_star))*(SFHdemo_nodust$flux[,2]-SFHdemo_dust$flux[,2]))

plot(SFHdemo_nodust$flux, log='xy', ylim=c(1e-20,1e-15), type='l',
xlab=BC03lr$Labels$Wavelab, ylab='Flux (erg/s/cm^2/Ang)')
lines(SFHdemo_dust$flux, col='grey')
lines(10^seq(5,7,by=0.01), greybody_norm(10^seq(5,7,by=0.01), z=0.1, norm=total_atten),
col='brown')
Dale_temp=Dale_interp(type='NormSFR')
lines(Dale_Msol$Wave, Dale_temp$Aspec*total_atten, col='red')

wave=sort(c(SFHdemo_dust$flux[,'wave'], Dale_Msol$Wave))
SpecAtten=approxfun(SFHdemo_dust$flux, rule=2)(wave)+
          approxfun(Dale_Msol$Wave, Dale_temp$Aspec, rule=2)(wave/1.1)
lines(wave, SpecAtten, col='darkgreen')

#Some different SFHs with identical stellar mass formed (10^10 Msol):

bursty=SFHfunc(m1=2, m2=1, m3=0, m4=0, m5=0, m6=0, forcemass=1e10)
constant=SFHfunc(m1=1, m2=1, m3=1, m4=1, m5=1, m6=1, forcemass=1e10)
young=SFHfunc(m1=0, m2=0, m3=1, m4=1, m5=0, m6=0, forcemass=1e10)
dead=SFHfunc(m1=0, m2=0, m3=0, m4=0, m5=1, m6=1, forcemass=1e10)

#SFHs:

plot(BC03lr$Age, bursty$SFR, log='xy', type='l', col='purple', ylim=c(1e-3,1e2), 
  xlab='Age / Yrs', ylab='SFR / Msol/yr')
lines(BC03lr$Age, constant$SFR, col='blue')
lines(BC03lr$Age, young$SFR, col='darkgreen')
lines(BC03lr$Age, dead$SFR, col='red')
legend('bottomleft', legend=c('bursty', 'constant', 'young', 'old'), lty=1, 
  col=c('purple', 'blue', 'darkgreen', 'red'))

#Luminosities with default dust attenuation:

wave=bursty$wave_lum
plot(wave, bursty$lum_atten, log='xy', type='l', col='purple', ylim=c(1e-1,1e8), 
  xlab='Wavelength / Ang', ylab='Lsun / Ang')
lines(wave, constant$lum_atten, col='blue')
lines(wave, young$lum_atten, col='darkgreen')
lines(wave, dead$lum_atten, col='red')
legend('topright', legend=c('bursty', 'constant', 'young', 'old'), lty=1, 
  col=c('purple', 'blue', 'darkgreen', 'red'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SFH }% use one of  RShowDoc("KEYWORDS")
