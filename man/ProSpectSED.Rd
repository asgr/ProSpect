\name{ProSpectSED}
\alias{ProSpectSED}
\alias{ProSpectSEDlike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Top Level ProSpect SED Functions
}
\description{
\code{ProSpectSED} computes a fully attenuated and re-emitted SED model incompassing star light (for a number of different star formation history models), dust and AGN. \code{ProSpectSEDlike} is an interface for computing likelihoods, and by extension is the workhorse code that allows model fitting using a number of techniques.
}
\usage{
ProSpectSED(SFH = SFHfunc, z = 0.1, tau_birth = 1, tau_screen = 0.3,
    tau_AGN = 1, pow_birth = -0.7, pow_screen = -0.7, pow_AGN = -0.7,
    alpha_SF_birth = 1, alpha_SF_screen = 3, alpha_SF_AGN = 0,
    AGNlum = 0, sparse = 5, speclib = NULL, Dale = NULL, AGN = NULL,
    filtout = NULL, filters = "all", Dale_M2L_func = NULL, returnall = TRUE,
    H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, waveout = seq(2,
        9.35, by = 0.01), ref, unimax = 1.38e+10, agemax = NULL,
    LumDist_Mpc = NULL, addradio_SF = FALSE, addradio_AGN = FALSE,
    Te_SF = 10000, ff_frac_SF = 0.1, ff_power_SF = -0.1, sy_power_SF = -0.8,
    Te_AGN = 10000, ff_frac_AGN = 0.1, ff_power_AGN = -0.1,
    sy_power_AGN = -0.8, AGNct = 60, AGNal = 4, AGNbe = -0.5,
    AGNta = 1, AGNrm = 60, AGNan = 30, Eb = 0, L0 = 2175.8, LFWHM = 470,
    IGMabsorb = 0, ...)

ProSpectSEDlike(parm = c(8,9,10,10,0,-0.5,0.2), Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SFH}{
Function; the star formation history to use. Must be one of \code{\link{SFHfunc}} of \code{\link{SFHburst}} (default).
}
  \item{z}{
Numeric scalar; cosmological redshift of the source. If the provided \option{z} = 0 then the returned fluxes will be absolute, i.e. absolute magnitudes if \option{outtype}='mag'. If \option{z} <= 0 then only luminosities are computed if \option{filtout}=NULL. If filters are provided via \option{filtout} then the \option{Photom} output returns Absolute AB magnitudes for convenience.
}
  \item{tau_birth}{
Numeric scalar; dust tau of birth clouds. Associated with sub 10 Myr star formation. See \code{\link{CF_birth}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF_screen}}.
}
  \item{tau_AGN}{
Numeric scalar; dust tau of the AGN. See \code{\link{CF}}.
}
  \item{pow_birth}{
Numeric scalar; power to further raise the dust attenuation of birth clouds. See \code{\link{CF_birth}}.
}
  \item{pow_screen}{
Numeric scalar; power to further raise the dust attenuation of the dust screen. See \code{\link{CF_screen}}.
}
  \item{pow_AGN}{
Numeric scalar; power to further raise the dust attenuation of the AGN torus. See \code{\link{CF_screen}}.
}
  \item{alpha_SF_birth}{
Numeric scalar; alpha slope of the birth cloud dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_screen}{
Numeric scalar; alpha slope of the screen dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_AGN}{
Numeric scalar; alpha slope of the AGN torus dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{AGNlum}{
Numeric scalar; AGN bolometric luminosity in erg/s.
}
  \item{sparse}{
Numeric scalar; amount of sparse sampling of the spectra to make.
}
  \item{speclib}{
Object; optional. Pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}. Default of NULL will load the \code{\link{BC03lr}} library, but faster to pass it in if being used for e.g. fitting.
}
  \item{Dale}{
Object or scalar logical; optional. Pass in the Dale dust library \code{\link{Dale_NormTot}} directly (otherwise will run, but not work correctly). Default of NULL will load the \code{\link{Dale_NormTot}} library, but faster to pass it in if being used for e.g. fitting. If FALSE it will not re-emit the attenuated light in the MIR/FIR using Dale templates. This might be desired if the range of interest is purely the UV-NIR and you want faster generation/fitting (since the re-emitted dust will not be notable until at least the MIR).
}
  \item{AGN}{
Object; optional. Pass in the AGN template to use. Must be one of \code{\link{AGN_UnOb}}, \code{\link{AGN_UnOb_Sparse}} or \code{\link{Fritz}}. The latter is preferred and uses the parametric Fritz model (see \code{\link{Fritz}}).
}
  \item{filtout}{
Object; required. Pass in the photometric filters directly (either matrices of functions). If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} or filtout=foreach(i = filters)\%do\%{approxfun(getfilt(i))}, where the latter will be faster (see \code{\link{getfilt}}).
}
  \item{filters}{
Character vector; names of filters to use. See \code{\link{filters}} for options. Default 'all' selects all filters, 'GAMA' returns just the classic GAMA survey filters. Can also be a list of matrices or functions to pass into \code{\link{bandpass}}. This creates increasingly faster code, but requires more user effort.
}
  \item{Dale_M2L_func}{
Function; if dust masses and luminosities are desired then the user should supply a function to convert \option{alpha_SF} to mass-to-light, e.g. \code{\link{Dale_M2L_func}}
}
  \item{returnall}{
Logical; if true then returns all the various sub SEDs. See Details. If FALSE then code just returns \option{photom_out} (photometry in Jansky).
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
  \item{waveout}{
Numeric vector; desired output log10 wavelength grid to use in Ang. Default covers typical range of galaxy SED data at a resolution good enough for broad band photometry. Spectroscopic work requires the resolution to be increase by a factor ~10. If \option{emission}=TRUE (to be passed down to \code{\link{SFHfunc}}) and \option{waveout} is missing (so not explicitly supplied by the user) then \option{waveout} is internally set to NULL, meaning the full available resolution is used to ensure a good reproduction of the emission features. This operation is generally desirable, since it would be a bad idea to produce emission features, but then sample the wavelength range too coarsely to see any!
}
  \item{ref}{
The name of a reference cosmology to use, one of 137 / 737 / Planck / Planck13 / Planck15 / Planck18 / WMAP / WMAP9 / WMAP7 / WMAP5 / WMAP3 / WMAP1 / Millennium / GiggleZ. Planck = Planck18 and WMAP = WMAP9. The usage is case insensitive, so wmap9 is an allowed input. This overrides any other settings for H0, OmegaM and OmegaL. If OmegaR is missing from the reference set then it is inherited from the function input (0 by default). See \code{\link{cosref}} for details.
}
  \item{unimax}{
Numeric scalar; maximum allowed age of any stellar population relative to z=0 (i.e. today) in years. Any star formation that is older than this will be set to 0. Overridden by \option{agemax} if that is provided.
}
  \item{agemax}{
Numeric scalar; maximum allowed age of any stellar population relative to the redshift of the observed galaxy in years. Any star formation that is older than this will be set to 0. Overridden by \option{unimax} if provided. It is computational a bit faster to provide this since we do not need to compute a lookup time from the redshift.
}
  \item{LumDist_Mpc}{
Numeric scalar; Luminosity distance computed in units of cm. Default is NULL. The luminosity distance can be supplied for repeated computations when the redshift of the object is constant. This should be supplied using the "correct" cosmology (as near as possible), i.e. using H0 ~ 70 km/s/Mpc (not 100).
}
  \item{addradio_SF}{
Logical scalar; should free-free and synchrotron radio continuum linked to star-formation be computed and appended to the spectrum? If TRUE, emission range is controlled by \option{waveout}. Uses the \code{\link{radiocont}} function to compute and append the radio continuum (see there for details). Note that even with this set to FALSE (default) the output spectrum will still have the radio continuum computed for the \code{\link{Dale}} templates. The default arguments almost perfectly reproduce the  \code{\link{Dale}} results, so it is only necessary to set \option{addradio_SF} = TRUE if you want to change the \option{Te}, \option{ff_frac}, \option{ff_power} and \option{sy_power} parameters to non-default values.
}
  \item{Te_SF}{
Numeric scalar; effective HII nebular plasma temperature in Kelvin used in \code{\link{radiocont}}. Leave as the default (10,000 K), unless you know what you are doing.
}
  \item{ff_frac_SF}{
Numeric scalar; free-free HII nebular plasma thermal radio emission contribution fraction at 1.4 GHz used in \code{\link{radiocont}} (should be between 0.05 and 0.2 typically). By definition the fraction of radio synchrotron emission associated with star-formation at 1.4 GHz will be 1 - \option{ff_frac_SF}.
}
  \item{ff_power_SF}{
Numeric scalar; power-law slope of the free-free nebular plasma thermal radio emission associated with star-formation used in \code{\link{radiocont}}. Leave as the default (-0.1), unless you know what you are doing.
}
  \item{sy_power_SF}{
Numeric scalar; power-law slope of the synchrotron radio emission associated with star-formation used in \code{\link{radiocont}} (should be between -0.8 and -0.6 typically).
}
  \item{addradio_AGN}{
Logical; should separate free-free and synchrotron radio continuum associated with an AGN be added? Uses the \code{\link{radiocont}} function to compute and append the radio continuum (see there for details). If FALSE (default) then any radio continuum computed for the AGN component within \code{\link{Dale}} will be subtracted regardless of the value of \option{addradio_SF}.
}
  \item{Te_AGN}{
Numeric scalar; effective HII nebular plasma temperature in Kelvin used in \code{\link{radiocont}} for the AGN component. Leave as the default (10,000 K), unless you know what you are doing.
}
  \item{ff_frac_AGN}{
Numeric scalar; free-free HII nebular plasma thermal radio emission contribution fraction at 1.4 GHz used in \code{\link{radiocont}} for the AGN component (should be between 0.05 and 0.2 typically). By definition the radio synchrotron emission fraction of the AGN component at 1.4 GHz will be 1 - \option{ff_frac_AGN}.
}
  \item{ff_power_AGN}{
Numeric scalar; power-law slope of the free-free nebular plasma thermal radio emission associated with an AGN used in \code{\link{radiocont}}. Leave as the default (-0.1), unless you know what you are doing.
}
  \item{sy_power_AGN}{
Numeric scalar; power-law slope of the synchrotron radio emission associated with an AGN used in \code{\link{radiocont}} (should be between -0.8 and -0.6 typically).
}
  \item{AGNct}{
Numeric scalar; opening angle of torus in degrees.
}
  \item{AGNal}{
Numeric scalar; gamma dust parameter.
}
  \item{AGNbe}{
Numeric scalar; Beta dust parameter.
}
  \item{AGNta}{
Numeric scalar; optical depth tau.
}
  \item{AGNrm}{
Numeric scalar; outer to inner torus radius ratio.
}
  \item{AGNan}{
Numeric scalar; angle of observation in degrees, where 0 is edge on with torus, and 90 is polar aligned.
}
  \item{parm}{
Numeric vector; arguments to be passed to \code{ProSpectSED}.
}
  \item{Data}{
Object; data object containing everything else needs. See Details.
}
  \item{Eb}{
Numeric scalar; Additional attenuation at the 2175.8 Ang dust bump. The default of 0 means this has no effect, MW is near 3.3, and typical galaxies are around \option{Eb}=1.
}
  \item{L0}{
Numeric scalar; location of the 2175.8 Ang dust bump in Angstroms (probably do not adjust this).
}
  \item{LFWHM}{
Numeric scalar; width of the 2175.8 Ang dust bump in Angstroms (probably do not adjust this).
}
  \item{IGMabsorb}{
Numeric scalar; absorbing fraction for the intervening IGM. This should be thought of as absorption due to material between the galaxy of interest an the observer. This suppressing fraction is applied to wavelength in the Lyman regime (911.8 - 1215.7 Angstrom). Based on Songaila (2004) a resonable approximation for this is given by \option{IGMabsorb} = pnorm(z, mean=3.8, sd=1.2), meaning it is almost 0 (complete transmission) at low redshift and almost 1 (complete absorption) by redshift 8.
}
  \item{\dots}{
All further arguments are passed directly into \option{SFH}, so by this mechanism you can set function specific arguments that control the star formation history etc.
}
}
\details{
This function allows full modelling of quite sophisticated star formation histories, dust attenuation and re-emission, as well as correctly attenuated AGN. It is useful directly, but it is also the workhorse function for the \code{\link{ProSpectSEDlike}} function, that allows us to compute likelihoods against target data and fit models in a generic manner.

Some more info on the required \option{Data} list inputs:

\describe{
  \item{Data$flux}{
Data.frame, required; 4 columns: filter [band name], cenwave [central wave, Ang], flux [Jy], fluxerr [Jy].
}
  \item{Data$SFH}{
Function, required; the star formation history function to use. Must be one of \code{\link{SFHfunc}} or \code{\link{SFHburst}}.
}
  \item{Data$arglist}{
List, optional; named arguments and values that are passed directly into \code{ProSpectSED}. Arguments that are not matched in \code{ProSpectSED} are passed directly onto \option{SFH}, so by this mechanism you can set function specific arguments that control the star formation history, i.e. this allows you to set values to non default even when not fitting these parameters. You almost certainly want to at least specify Data$arglist$z (the object redshift) since this almost always differ from the fiducial default of \option{z}=0.1.
}
  \item{Data$speclib}{
List, optional but recommended; pass in the spectral library directly. Must be one of \code{\link{BC03lr}}, \code{\link{BC03hr}}, \code{\link{EMILES}}.
}
  \item{Data$Dale}{
List, optional but recommended; pass in the Dale dust library directly (should be one of \option{Dale_NormTot}, \option{Dale_NormSFR} otherwise will run, but not work correctly). If set to FALSE then there will be no FIR re-emission.
}
  \item{Data$AGN}{
Data.frame, optional; pass in the AGN template to use. Must be one of \code{\link{AGN_UnOb}} or \code{\link{AGN_UnOb_Sparse}}.
}
  \item{Data$filtout}{
List; required. Pass in the photometric filters directly. If \code{filters} is a vector of filter names, this can be achieved efficiently with a command like: filtout=foreach(i = filters)\%do\%{getfilt(i)} or filtout=foreach(i = filters)\%do\%{approxfun(getfilt(i))}, where the latter will be faster (see \code{\link{getfilt}}).
}
  \item{Data$fit}{
Character scalar, required; one of 'optim', 'LD' or 'check'. See Value for full information on what is output in each case. If missing it is set to 'optim'.
}
  \item{Data$like}{
Character scalar, required; specifies the type of likelihood. One of 'norm' (Normal), 'st' (Student-T) or 'chisq' (Chi-Square). If not one of these then all useful output is returned in a list. If missing it is set to 'st'.
}
  \item{Data$parm.names}{
Character vector, required; correctly named and position matched parameters for fitting. Must match length of \option{parm} in \code{ProSpectSEDlike}, and elements should correctly describe the respective values.
}
  \item{Data$logged}{
Logical vector, optional; match length of \option{parm} in \code{ProSpectSEDlike}, and elements should say which values should be raised to the power 10 when fitting. This is useful since scale parameters are better fit in log space, so Data$logged=TRUE for those elements. If missing assumes all parameters are to be fit in linear space.
}
  \item{Data$prior}{
Function, optional; a prior function that takes in the \option{parm} vector and returns a log-likelihood to be added to LL to form LP. If not provided then it is ignored internally. The user can expect the vector \option{parm} to have elements specified by \option{Data$parm.names} for ease.
}
  \item{Data$intervals}{
List, optional; has elements \option{lo} and \option{hi} which if provided must be the same length as \option{parm} in \code{ProSpectSEDlike}. These specify the lower and upper limits in the same linear/log space as \option{parm}, i.e. set \option{lo} to -1 if you do not want to logged unit to go below 0.1.
}
  \item{Data$constraints}{
Function, optional; a constraints function that takes in the\option{parm} vector and modifies it as required, updating with the new \option{parm} which is passed out. Using this you can control physical relationships between certain parameters if required. \option{Data$intervals} is passed after this, and always over-ride the constraints if there is a conflict. The user can expect the vector \option{parm} to have elements specified by \option{Data$parm.names} for ease.
}
  \item{Data$mon.names}{
Character vector, optional; "LP" (log posterior), "masstot" (total stellar mass formed), "massrem" (total stellar mass remaining) and "SFRburst" (stars formed in the most recent 1e8 years). For dust can contain "dustmass.birth", "dustmass.screen", "dustmass.AGN", "dustmass.total" (dust masses in birth cloud, screen, AGN and total), "dustlum.birth", "dustlum.screen", "dustlum.AGN", "dustlum.total" (dust luminosities in birth cloud, screen, AGN and total). You can also monitor fluxes via "flux.XXX", where XXX specifies the name of the band as given in Data$flux$filter. Only named monitored variables are returned, in the order specified.
}
  \item{Data$logbin}{
Logical, optional; should the automatic wavelength grid be evenly log spaced? If TRUE it will be evenly log spaced, if FALSE it will be evenly linearly spaced. Passed into \code{\link{specReBin}}. Say you have a wavelength 10 and 20, a linear bin would be at 15 and a log bin would be at 14.1.
}
  \item{Data$rough}{
Logical, otional; should rough interpolation be used to rebin the spectrum (default is TRUE) or should full be pixel integration be done (FALSE). Passed into \code{\link{specReBin}}.
}
  \item{Data$verbose}{
Logical, optional; should verbose messages be shown. If missing it is set to TRUE.
}
}
}
\value{
For \code{ProSpectSED}: \cr

If \option{returnall}=FALSE then we get a vector containing the \option{Photom} output, which is usually the photometry in Jansky (special case when \option{z} <= 0, where this is instead the Absolute AB magnitudes for convenience).

If \option{returnall}=TRUE, then we get a list containing:

\item{Photom}{Vector; photometry in Jansky. If \option{z} <= 0 and filters are provided via \option{filtout} then the \option{Photom} output returns Absolute AB magnitudes for convenience (i.e. source is placed at a distance of 10 pc).}
\item{FinalFlux}{Two column data.frame; Wavelength in Ang (wave), Flux density in Jansky (flux)}
\item{FinalLum}{Two column data.frame; Wavelength in Ang (wave), Total SED luminosity in Lsol / Ang (lum)}
\item{StarsAtten}{Two column data.frame; Wavelength in Ang (wave), Attenuated stellar SED luminosity in Lsol / Ang (lum)}
\item{StarsUnAtten}{Two column data.frame; Wavelength in Ang (wave), Unattenuated stellar SED luminosity in Lsol / Ang (lum)}
\item{DustEmit}{Two column data.frame; Wavelength in Ang (wave), Re-emitted stellar dust SED luminosity in Lsol / Ang (lum)}
\item{AGN}{Two column data.frame; Wavelength in Ang (wave), Attenuated and re-emitted AGN SED luminosity in Lsol / Ang (lum)}
\item{Stars}{List; the full output of the relevant SFH function used. See \code{\link{SFHfunc}}}
\item{dustmass}{Numeric vector; contains dust mass for the birth cloud, screen, AGN and total in Msol (requires \option{Dale_M2L_func})}
\item{dustlum}{Numeric vector; contains dust luminosities for the birth cloud, screen, AGN and total in Lsol (requires \option{Dale_M2L_func})}
\item{call}{The original function call.}
\item{z}{Numeric scalar; the input redshift.}
\item{filtout}{List; the \option{filtout} list used.}
\item{filters}{List; the \option{filters} list used.}
\item{cosmo}{List; the input 'H0', 'OmegaM' and 'OmegaL'.}

For \code{ProSpectSEDlike}: \cr

Output depends on the flag for \option{Data$fit}, where it can be one of 'optim' (or 'CMA'), 'LD' (or 'LA'), or 'check'. \cr

The 'optim' output returns just the log-posterior value (LP). \cr

The 'LD' output returns the full structure required for \code{LaplacesDemon} and \code{LaplaceApproximation} and should be set as such when using those optimisers. \cr

The 'check' output return all the outputs for \code{LaplacesDemon} and \code{LaplaceApproximation} but additionally the full output of \code{ProSpectSED} with \option{returnall}=TRUE set as list item 'SEDout', and the input \option{Data} object as list item 'Data'. If you have the \code{ParmOff} package (on Github asgr/ParmOff, required for correct argument matching between the various functions used) it will also output the various stellar mass calculations of \code{\link{SMstarfunc}} including stellar mass remaining (TotSMstar) which is the reason this function is usually run.

In the last case of \option{Data$fit}='check' it will be returned with class 'ProSpectSEDlike', which allows you to create a quick diagnostic plot via the \code{plot} method.
}

\references{
Robotham et al, 2020, MNRAS, 495, 905 \cr
Songaila, 2004, AJ, 127, 2598
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{SFH}}, \code{\link{magABcalc}}, \code{\link{dustmass}}, \code{\link{optim}}
}

\examples{
data('BC03lr')
data('Dale_NormTot')
data('AGN_UnOb_Sparse')
data('Dale_M2L_func')
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
filtout={}
for(i in filters){filtout=c(filtout,list(approxfun(getfilt(i))))}
#for(i in filters){filtout=c(filtout,list(getfilt(i)))} #This works too, but slower!

#Try playing around with the below to see what impact different AGN, dust etc have:

testSED=ProSpectSED(AGNlum=1e43, tau_birth=1, tau_screen=0.3, tau_AGN=2,
alpha_SF_birth=1, alpha_SF_screen=3, alpha_SF_AGN=0, speclib=BC03lr, Dale=Dale_NormTot,
AGN=AGN_UnOb_Sparse, filtout=filtout, Dale_M2L_func=Dale_M2L_func, returnall=TRUE)

plot(testSED$FinalLum, log='xy', xlim=c(1e2,1e7), ylim=c(1e2,1e7),
  xlab='Wavelength (Ang)', ylab='Lum (Lsol/Ang)', type='l', lwd=5)
lines(testSED$StarsUnAtten, col='blue', lty=2)
lines(testSED$StarsAtten, col='green')
lines(testSED$DustEmit, col='brown')
lines(testSED$AGN, col='purple')
legend('topright',
legend=c('Total Lum', 'Star Un-Atten', 'Stars Atten', 'Dust Emit', 'AGN'),
  col=c('black', 'blue', 'darkgreen', 'brown', 'purple'),
  lty=c(1,2,1,1,1),
  lwd=c(5,1,1,1,1)
)

# or easier:

plot(testSED)

#An example of a more physical SFH and ZH:

testSED2 = ProSpectSED(massfunc = massfunc_snorm, Z=Zfunc_massmap_box)
plot(testSED2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SED }% use one of  RShowDoc("KEYWORDS")
