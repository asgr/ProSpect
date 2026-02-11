\name{Dale_M2L_variableDTH_func}
\alias{Dale_M2L_variableDTH_func}
\title{
Dale Template Mass-to-Light variable DTH
}
\description{
Returns the correct mass-to-light for a given \option{alpha_SF} (which is akin to dust temperature) similar to \link[ProSpectData]{Dale_M2L}.
}
\usage{
Dale_M2L_variableDTH_func(alpha_SF, qPAH_VSG = 0.14, pivwave = 10^5.4, step_speed = 9)
}
\arguments{
  \item{alpha_SF}{Numeric vector; desired interpolated alpha slope of the star forming population. Lower values mean hotter dust.}
  \item{qPAH_VSG}{Numeric; between 0 and 1. The mass fraction in very-small-grains and PAH molecules, recommended 0.1-0.14. }
  \item{pivwave}{Numeric; angstroms. Wavelength at which the dust emission transitions from very-small-grains+PAHs (IR+MIR) to big-grains (FIR).}
  \item{step_speed}{Numeric. Controls the sharpness of the step-wise transition from very-small-grains+PAHs to big-grains.}
}
\details{
This function allows flexibility for the fraction of the dust mass contributed to by both very-small-grains and PAH molecules (which emit in the IR-MIR) compared to big grains (which emit in the FIR). For further details, see Fig. 4 and accompyaning discussion in D'Silva+26 (https://arxiv.org/abs/2601.08112).
}
\value{
Numeric vector; M2L outputs as a function of \option{alpha_SF}.
}
\references{
D'Silva+26 (https://arxiv.org/abs/2601.08112)
}
\author{
Jordan D'Silva
}
\examples{
plot(
  seq(0.06, 4, 0.1),
  Dale_M2L_variableDTH_func(seq(0.06, 4, 0.1)),
  log='y',
  xlim = c(0.05, 4.05),
  ylim = c(50, 1e6),
  xlab = "alpha_SF",
  ylab = "Dale M2L"
)
points(
  seq(0.06, 4, 0.1),
  Dale_M2L_variableDTH_func(seq(0.06, 4, 0.1), qPAH_VSG = 0.0, pivwave = 1),
  pch = 16
)
legend(
  x = "topright",
  pch = c(1,16),
  legend = c("Variable DTH", "qPAH_VSG=0.0, pivwave=1 \n -> same as Dale_M2L")
)

prospect_old_dust_mass = ProSpectSED(
  Dale = Dale_NormTot,
  Dale_M2L_func = Dale_M2L_func
)
prospect_new_dust_mass = ProSpectSED(
  Dale = Dale_NormTot,
  Dale_M2L_func = Dale_M2L_variableDTH_func
)

print(
  prospect_old_dust_mass$dustmass
)
print(
  prospect_new_dust_mass$dustmass
)
print(
  prospect_old_dust_mass$dustmass/prospect_new_dust_mass$dustmass
) # factor of ~3 lower dust mass...

prospect_old_dust_mass_again = ProSpectSED(
  Dale = Dale_NormTot,
  Dale_M2L_func = function(alpha_SF){
    Dale_M2L_variableDTH_func(
      alpha_SF,
      qPAH_VSG = 0.0,
      pivwave = 1
    )
  }
)

print(
  prospect_old_dust_mass_again$dustmass
) ## same as Dale_M2L
}
