# Sanity checks that the model responds to physical parameters in the direction
# physics demands. These are the primary tests in the suite: they are
# deliberately qualitative (directions and monotonicity) so that they catch
# wiring bugs without pinning down the exact numerics of any SSP library.

test_that('dust attenuation reddens the SED monotonically in optical/NIR bands', {
  taus = c(0, 0.25, 0.5, 1, 2, 4)

  photom = vapply(taus, function(tau) {
    .ps_photom(tau_screen = tau, tau_birth = 0)
  }, numeric(length(.ps_bands)))
  rownames(photom) = .ps_bands

  # More screen dust always removes more starlight from the bands that are
  # dominated by starlight. W1/W2 are excluded: the Dale dust continuum is
  # already bright there, so they can go either way.
  for (band in .ps_stellar_bands[1:9]) {
    expect_true(
      all(diff(photom[band, ]) < 0),
      label = paste0(band, ' decreases with tau_screen')
    )
  }

  # ...and it removes far more FUV than NIR light, so the SED reddens.
  attenuation_fuv = photom['FUV_GALEX', 6] / photom['FUV_GALEX', 1]
  attenuation_k = photom['K_VISTA', 6] / photom['K_VISTA', 1]
  expect_lt(attenuation_fuv, 1e-2)
  expect_gt(attenuation_k, 0.1)
  expect_lt(attenuation_fuv, attenuation_k)
})

test_that('dust optical depth raises FIR emission (energy is reradiated)', {
  # dustlum is only populated when the dust mass-to-light function is supplied,
  # so use the emitted spectrum itself rather than the summary field.
  taus = c(0, 0.25, 0.5, 1, 2, 4)

  for (param in c('tau_screen', 'tau_birth')) {
    emitted = vapply(taus, function(tau) {
      args = list(tau_birth = 0, tau_screen = 0)
      args[[param]] = tau
      .ps_integral(do.call(.ps_sed, args)$DustEmit)
    }, numeric(1))

    expect_true(
      all(diff(emitted) > 0),
      label = paste0('dust emission increases with ', param)
    )
    expect_equal(emitted[1], 0, tolerance = 1e-12)
  }

  # The FIR/sub-mm bands are where this shows up in observable photometry.
  fir_growth = vapply(taus, function(tau) {
    unname(.ps_photom(tau_screen = tau, tau_birth = 0)['S250_Herschel'])
  }, numeric(1))
  expect_true(all(diff(fir_growth) > 0))
})

test_that('birth-cloud dust is more damaging to UV than screen dust', {
  # A birth cloud only enshrouds stars younger than ~10 Myr, so the same tau
  # should suppress the FUV far less than an equal tau in the screen (which
  # attenuates every stellar population), while leaving the NIR nearly alone.
  clean = .ps_photom(tau_birth = 0, tau_screen = 0)
  cloudy_birth = .ps_photom(tau_birth = 4, tau_screen = 0)
  cloudy_screen = .ps_photom(tau_birth = 0, tau_screen = 4)

  fuv_birth = unname(cloudy_birth['FUV_GALEX'] / clean['FUV_GALEX'])
  fuv_screen = unname(cloudy_screen['FUV_GALEX'] / clean['FUV_GALEX'])
  k_birth = unname(cloudy_birth['K_VISTA'] / clean['K_VISTA'])

  expect_gt(fuv_birth, fuv_screen)
  expect_gt(fuv_birth, 0.1)   # birth cloud alone leaves a UV floor...
  expect_lt(fuv_screen, 1e-3) # ...the screen does not
  expect_gt(k_birth, 0.9)     # and neither touches the NIR much
})

test_that('dust absorption is exactly balanced by dust re-emission', {
  # The Dale templates are normalised to unit total luminosity, so absorbed
  # starlight must equal the emitted dust luminosity. This is a strong
  # invariant on the attenuation + emission coupling.
  seds = list(
    .ps_sed(tau_birth = 0.5, tau_screen = 0.5),
    .ps_sed(tau_birth = 2, tau_screen = 0.2),
    .ps_sed(tau_birth = 0.1, tau_screen = 3)
  )

  for (sed in seds) {
    absorbed = sed$Stars$lumtot_birth + sed$Stars$lumtot_screen
    emitted = .ps_integral(sed$DustEmit)

    expect_gt(absorbed, 0)
    # Note lumtot_atten in the Stars list is the *absorbed* luminosity, not the
    # transmitted one; it equals lumtot_birth + lumtot_screen by construction.
    expect_equal(sed$Stars$lumtot_atten, absorbed, tolerance = 1e-10)
    expect_equal(emitted, absorbed, tolerance = 1e-8)

    # And transmitted starlight is the unattenuated spectrum minus the absorbed part.
    transmitted = .ps_integral(sed$StarsAtten)
    expect_equal(transmitted, sed$Stars$lumtot_unatten - absorbed, tolerance = 1e-8)
  }
})

test_that('a dust-free model leaves the stellar spectrum unattenuated', {
  sed = .ps_sed(tau_birth = 0, tau_screen = 0)

  # lumtot_atten is the absorbed luminosity, so with no dust it must vanish.
  expect_equal(sed$Stars$lumtot_atten, 0, tolerance = 1e-12)
  expect_equal(sed$Stars$lumtot_birth, 0, tolerance = 1e-12)
  expect_equal(sed$Stars$lumtot_screen, 0, tolerance = 1e-12)
  expect_true(all(sed$DustEmit$lum == 0))
  expect_equal(.ps_integral(sed$DustEmit), 0, tolerance = 1e-12)
  expect_equal(sed$Stars$lum_atten, sed$Stars$lum_unatten, tolerance = 1e-12)
})

test_that('higher metallicity reddens and dims the UV', {
  # Cooler, more metal-rich stars are redder: the UV-to-NIR ratio must fall as
  # Z rises across the metallicities BC03 actually provides.
  Zs = c(0.0001, 0.004, 0.02, 0.05)

  fuv_k = vapply(Zs, function(Z) {
    .ps_ratio(.ps_photom(Z = Z), 'FUV_GALEX', 'K_VISTA')
  }, numeric(1))
  g_k = vapply(Zs, function(Z) {
    .ps_ratio(.ps_photom(Z = Z), 'g_SDSS', 'K_VISTA')
  }, numeric(1))

  expect_true(all(diff(fuv_k) < 0))
  expect_true(all(diff(g_k) < 0))
})

test_that('a metallicity history interpolates between its endpoints', {
  # Zfunc_massmap_lin maps the cumulative mass to [Zstart, Zfinal]; raising
  # Zfinal must redden the SED and raise the delivered Zvec.
  zfinal = c(0.001, 0.02, 0.04)

  sed = lapply(zfinal, function(Zf) {
    .ps_sed(Z = Zfunc_massmap_lin, Zfinal = Zf)
  })
  fuv_k = vapply(zfinal, function(Zf) {
    .ps_ratio(.ps_photom(Z = Zfunc_massmap_lin, Zfinal = Zf), 'FUV_GALEX', 'K_VISTA')
  }, numeric(1))

  expect_true(all(diff(fuv_k) < 0))

  for (i in seq_along(zfinal)) {
    expect_gte(min(sed[[i]]$Stars$Zvec), 1e-4)
    expect_lte(max(sed[[i]]$Stars$Zvec), zfinal[i] + 1e-6)
  }
  expect_lt(max(sed[[1]]$Stars$Zvec), max(sed[[3]]$Stars$Zvec))
})

test_that('a young burst makes the SED bluer and more luminous in the UV', {
  # massfunc_b5 phase 1 is 0-100 Myr. Holding the older phases fixed at 1, a
  # larger young phase adds hot stars, so FUV/K rises monotonically.
  m1 = c(0, 0.5, 3, 30)

  fuv_k = vapply(m1, function(m) {
    .ps_ratio(.ps_photom(m1 = m, m2 = 1), 'FUV_GALEX', 'K_VISTA')
  }, numeric(1))

  expect_true(all(diff(fuv_k) > 0))
  expect_gt(fuv_k[4] / fuv_k[1], 10)
})

test_that('SED flux scales linearly with the total mass formed', {
  # Nothing in the pipeline is non-linear in mass, so tripling m1 (which scales
  # every phase of massfunc_b5) must triple every photometric point.
  base = .ps_photom(m1 = 1)
  tripled = .ps_photom(m1 = 3)

  expect_equal(unname(tripled / base), rep(3, length(base)), tolerance = 1e-10)

  # and masstot must follow, while the mass-to-light ratio must not change.
  expect_equal(.ps_sed(m1 = 3)$Stars$masstot, 3 * .ps_sed(m1 = 1)$Stars$masstot,
               tolerance = 1e-10)
})

test_that('redshift dims and reddens the observed SED', {
  zs = c(0.01, 0.05, 0.1, 0.3, 0.5)

  photom = vapply(zs, function(z) .ps_photom(z = z), numeric(length(.ps_bands)))
  rownames(photom) = .ps_bands

  # Surface brightness dimming wins over the (1+z) band shift: flux falls.
  expect_true(all(diff(photom['K_VISTA', ]) < 0))

  # K-correction reddens the optical/NIR colour up to z = 0.5.
  g_k = photom['g_SDSS', ] / photom['K_VISTA', ]
  expect_true(all(diff(g_k) < 0))
})

test_that('cosmology changes the inferred flux at fixed redshift', {
  # Flux scales as 1/D_L^2, so a lower H0 (larger D_L) must give a fainter SED.
  H0s = c(50, 67.8, 100)

  k_band = vapply(H0s, function(H0) unname(.ps_photom(z = 0.3, H0 = H0)['K_VISTA']),
                  numeric(1))

  expect_true(all(diff(k_band) > 0))
})

test_that('LumDist_Mpc overrides the internal cosmology exactly', {
  z = 0.3
  DL = celestial::cosdistLumDist(z = z, H0 = 67.8, OmegaM = 0.308)

  computed = .ps_photom(z = z)
  supplied = .ps_photom(z = z, LumDist_Mpc = DL)
  expect_equal(unname(supplied / computed), rep(1, length(computed)), tolerance = 1e-10)

  # Halving the luminosity distance must quadruple the flux (1/D_L^2).
  closer = .ps_photom(z = z, LumDist_Mpc = DL / 2)
  expect_equal(unname(closer / computed), rep(4, length(computed)), tolerance = 1e-9)
})

test_that('Dale radiation-field strength softens the FIR SED', {
  # Higher alpha_SF means a harder radiation field and hence warmer dust, which
  # shifts the FIR peak blueward. P100/S250 falls monotonically.
  alphas = c(0.0625, 0.25, 0.5, 1, 2, 4)

  p100_s250 = vapply(alphas, function(a) {
    .ps_ratio(.ps_photom(alpha_SF_screen = a), 'P100_Herschel', 'S250_Herschel')
  }, numeric(1))

  expect_true(all(diff(p100_s250) < 0))
})

test_that('AGN luminosity raises mid-IR and UV flux monotonically', {
  e = new.env()
  data('Fritz', envir = e)

  lums = c(1e42, 1e43, 1e44)

  w3 = vapply(lums, function(lum) {
    unname(.ps_photom(AGN = e$Fritz, AGNlum = lum)['W3_WISE'])
  }, numeric(1))
  fuv = vapply(lums, function(lum) {
    unname(.ps_photom(AGN = e$Fritz, AGNlum = lum)['FUV_GALEX'])
  }, numeric(1))

  expect_true(all(diff(w3) > 0))
  expect_true(all(diff(fuv) > 0))
})

test_that('IGM absorption suppresses rest-frame FUV flux at high redshift', {
  clean = .ps_photom(z = 3, IGMabsorb = 0)
  absorbed = .ps_photom(z = 3, IGMabsorb = 'Inoue14')

  expect_lt(unname(absorbed['FUV_GALEX']), unname(clean['FUV_GALEX']))

  # The effect is a rest-frame UV cut-off, so the reddest band is untouched.
  expect_equal(unname(absorbed['K_VISTA']), unname(clean['K_VISTA']),
               tolerance = 1e-10)
})

test_that('SFHburst scales linearly with burst mass and keeps its colour', {
  masses = c(1e7, 1e8, 1e9)

  photom = vapply(masses, function(bm) {
    out = ProSpectSED(SFH = SFHburst, burstmass = bm, speclib = .ps_speclib(),
                      Dale = .ps_dale(), filtout = .ps_filtout(), filters = NULL,
                      z = 0.1, sparse = 10)
    stats::setNames(out$Photom, .ps_bands)
  }, numeric(length(.ps_bands)))
  rownames(photom) = .ps_bands

  # Colour is independent of normalisation...
  expect_equal(unname(photom[, 2] / photom[, 1]), rep(10, length(.ps_bands)),
               tolerance = 1e-10)
  # ...and a pure burst is much bluer than the default extended SFH.
  expect_gt(photom['FUV_GALEX', 1] / photom['K_VISTA', 1],
            .ps_ratio(.ps_photom(), 'FUV_GALEX', 'K_VISTA'))
})