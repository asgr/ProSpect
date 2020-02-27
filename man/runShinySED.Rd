\name{runShiny}
\alias{runShiny}
\alias{runShinySED}
\alias{runShinyFilters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Runs Shiny SED App
}
\description{
Simple wrapper function to run the Shiny SED and Filters App. The SED app shows model flux and luminosity. In the former case it is possible to overlay user data using the \option{flux} arguement. In the latter case you can see the individual components that add together to make the final SED. Note flux is per unit Hz (as standard for e.g. Jansky etc) and luminosity is per unit Angstrom (as standard for spectral libraries). The Filters app lets you convert between photometric filters using \code{\link{filterTranMags}}.
}
\usage{
runShinySED(flux = NULL, burstmass = 1e+08, youngmass = 1e+09, midmass = 1e+10,
oldmass = 1e+10, ancientmass = 1e+10, z = 0.1, tau_birth = 1, tau_screen = 0.3,
tau_AGN = 1, alpha_SF_birth = 1, alpha_SF_screen = 3, alpha_SF_AGN = 0, AGNlum = 1e+42,
Z = rep(5,5))
runShinyFilters()
}
\details{
Mostly used to run the app, but it can also return the current parameters to help with finding starting parameters for fitting.
}
\arguments{
  \item{flux}{
Data.frame(4 columns, filter name, central wave (Ang), flux (Jy), fluxerr (Jy)). If provided this data will be overlaid in the "Flux" tab to compare against the current model. See \code{\link{ProSpectSEDlike}}.
}
  \item{burstmass}{
Numeric scalar; the total stellar mass formed in a recent burst in units Msol. See \code{\link{SFHp5}}.
}
  \item{youngmass}{
Numeric scalar; the total stellar mass formed in young stars in units Msol. See \code{\link{SFHp5}}.
}
  \item{midmass}{
Numeric scalar; the total stellar mass formed in mid stars in units Msol. See \code{\link{SFHp5}}.
}
  \item{oldmass}{
Numeric scalar; the total stellar mass formed in old stars in units Msol. See \code{\link{SFHp5}}.
}
  \item{ancientmass}{
Numeric scalar; the total stellar mass formed in ancient stars in units Msol. See \code{\link{SFHp5}}.
}
  \item{z}{
Numeric scalar; cosmological redshift of the source. See \code{\link{SFHp5}}.
}
  \item{tau_birth}{
Numeric scalar; dust tau of birth clouds. Associated with sub 10 Myr star formation. See \code{\link{CF_birth}}.
}
  \item{tau_screen}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF_screen}}.
}
  \item{tau_AGN}{
Numeric scalar; dust tau of the dust screen. See \code{\link{CF}}.
}
  \item{alpha_SF_birth}{
Numeric scalar; alpha slope of the birth cloud dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_screen}{
Numeric scalar; alpha slope of the screen dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{alpha_SF_AGN}{
Numeric scalar; alpha slope of the AGN taurus dust. Lower values mean hotter dust. See \code{\link{Dale}}.
}
  \item{AGNlum}{
Numeric scalar; AGN bolometric luminosity in erg.
}
  \item{Z}{
Integer vector; length 5. Specifies the metallicity to use for each phase (default is 5, so solar metallicity for all phases, which is fairly appropriate at low \option{z}. See \code{\link{SFHp5}}.
}
}
\value{
\code{runShinySED} runs a Shiny app run, which in turn creates an interactive SED plot. The first tab shows the output in flux space for the full SED (all components combined) and the second tab shows the output in luminosity space (split by component also). The luminosity plane also shows the effective star formation history. The redded out region shows the ignored region that is older than the age of the Universe (13.8 Gyrs).

If the user hits \option{Stop!} in the app, then the current parameters are returned out in a list. This is useful for interactively finding good starting arguments for fitting purposes. The logged parameters in the Shiny app (masses and AGN luminosity) are returned in linear units since this means they can be directly used in \code{\link{ProSpectSED}}.

\code{runShinyFilters} lets you interactively convert between target and reference filters.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{ProSpectSED}}, \code{\link{filterTranMags}}
}
\examples{
\dontrun{
runShinySED()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ Shiny }% use one of  RShowDoc("KEYWORDS")
\concept{ SED }% __ONLY ONE__ keyword per line
\concept{ Filters }
