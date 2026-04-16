library(ProSpect)

wave <- seq(1200, 10000, length.out = 200)
tau <- 0.6

# Backward compatibility: CF path must be identical.
cf_old <- CF_screen(wave = wave, tau = tau, pow = -0.7, Eb = 0.8, L0 = 2175.8, LFWHM = 470)
cf_new <- screen_atten(
  wave = wave,
  tau = tau,
  pow = -0.7,
  Eb = 0.8,
  L0 = 2175.8,
  LFWHM = 470,
  dust_law = "CF"
)
stopifnot(isTRUE(all.equal(cf_old, cf_new, tolerance = 0)))

# Salim+18 basic behavior and pivot normalization.
salim <- Salim18_screen(wave = wave, tau = tau, delta = 0.2, B = 1)
pivot_val <- Salim18_screen(wave = 5500, tau = tau, delta = 0.2, B = 1)
stopifnot(abs(pivot_val - exp(-tau)) < 1e-12)
stopifnot(all(salim > 0 & salim <= 1))

# Eb fallback maps to B when using Salim18 via wrapper.
salim_from_Eb <- screen_atten(wave = wave, tau = tau, dust_law = "Salim18", delta = -0.3, Eb = 0.7)
salim_from_B <- screen_atten(wave = wave, tau = tau, dust_law = "Salim18", delta = -0.3, B = 0.7)
stopifnot(isTRUE(all.equal(salim_from_Eb, salim_from_B, tolerance = 1e-12)))

flux <- rep(1, length(wave))
atten <- CF_screen_atten(wave = wave, flux = flux, tau = tau, dust_law = "Salim18", delta = 0.1, B = 0.5)
stopifnot(atten$attenfrac >= 0 && atten$attenfrac <= 1)

cat("All Salim+18 dust attenuation tests passed.\n")
