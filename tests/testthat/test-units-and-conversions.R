# Unit-handling tests. ProSpect mixes Angstrom and Hz flux densities, Lsun/Ang
# and Jansky, and log10 and linear wavelengths. These tests pin down each
# convention so a units slip is caught here rather than in a fit.

test_that('AB magnitude and flux-density conversions are mutually inverse', {
  mags = c(-10, 0, 15, 20, 25, 30)

  expect_equal(Jansky2magAB(magAB2Jansky(mags)), mags, tolerance = 1e-12)
  expect_equal(CGS2magAB(magAB2CGS(mags)), mags, tolerance = 1e-12)
  expect_equal(magAB2Jansky(mags) * 1e-23, magAB2CGS(mags), tolerance = 1e-12)

  # The defining AB zeropoint: 3631 Jy at AB = 0.
  expect_equal(magAB2Jansky(0), 3630.78, tolerance = 1e-5)
  expect_equal(Jansky2magAB(1), 8.9, tolerance = 1e-12)
})

test_that('wavelength and frequency flux densities convert consistently', {
  # convert_wave2freq / convert_freq2wave are exact inverses on a spectrum.
  wave = seq(3000, 10000, by = 25)
  f_wave = exp(-((log10(wave) - 3.6)^2) / 0.02)

  f_freq = convert_wave2freq(f_wave, wave)
  expect_equal(convert_freq2wave(f_freq, wave), f_wave, tolerance = 1e-12)

  # The conversion is f_nu = f_lambda * lambda^2 / c, so the ratio is checkable
  # in closed form.
  expect_equal(f_freq / f_wave, wave^2 / .c_to_mps * 1e-10, tolerance = 1e-10)
})

test_that('Janskycalc and magABcalc reproduce the AB definition of a flat f_nu', {
  # A flat 1e-14 erg/s/cm^2/Hz spectrum is by construction 1e9 Jy.
  wave = seq(4000, 7000, by = 10)
  f_nu = rep(1e-14, length(wave))
  f_wave = convert_freq2wave(f_nu, wave)

  expect_equal(Janskycalc(wave, f_wave, 'r_SDSS'), 1e-14 * 1e23, tolerance = 1e-10)

  # The cgs AB zeropoint is -48.6 (f_nu in erg/s/cm^2/Hz), not the +8.9 used
  # for the Jansky form.
  expect_equal(magABcalc(wave, f_wave, 'r_SDSS'),
               -2.5 * log10(1e-14) - 48.6, tolerance = 1e-10)
})

test_that('Lum2Flux and Flux2Lum are inverses and redshift the wavelength', {
  wave = 10^seq(2.5, 4.5, by = 0.01)
  lum = exp(-((log10(wave) - 3.5)^2) / 0.02) * 1e10

  flux = Lum2Flux(wave, lum, z = 0.3, H0 = 67.8, OmegaM = 0.308)
  expect_equal(flux$wave, wave * 1.3, tolerance = 1e-10)

  back = Flux2Lum(flux$wave, flux$flux, z = 0.3, H0 = 67.8, OmegaM = 0.308)
  expect_equal(back$lum, lum, tolerance = 1e-10)
  expect_equal(back$wave, wave, tolerance = 1e-10)
})

test_that('cosmological dimming follows the inverse square law exactly', {
  wave = seq(4000, 8000, by = 100)
  lum = rep(1e10, length(wave))

  DL = celestial::cosdistLumDist(z = 0.3, H0 = 67.8, OmegaM = 0.308)

  # Supplying the distance directly must match computing it internally.
  from_ref = Lum2Flux(wave, lum, z = 0.3, H0 = 67.8, OmegaM = 0.308)
  from_supplied = Lum2Flux(wave, lum, z = 0.3, LumDist_Mpc = DL)
  expect_equal(from_supplied$flux, from_ref$flux, tolerance = 1e-10)

  # Doubling the distance quarters the flux.
  doubled = Lum2Flux(wave, lum, z = 0.3, LumDist_Mpc = 2 * DL)
  expect_equal(doubled$flux, from_ref$flux / 4, tolerance = 1e-10)
})

test_that('ProSpectSED Photom is Jansky at z > 0 and AB mag at z <= 0', {
  bright = .ps_photom(z = 0.1)
  absolute = .ps_photom(z = 0)

  # Jansky values at z = 0.1 for a normal galaxy are tiny; AB mags are O(10-25).
  expect_true(all(bright > 0 & bright < 1))
  expect_true(all(absolute > -30 & absolute < 0))

  # The two are related by the distance modulus between z = 0 and z = 0.1.
  expect_false(isTRUE(all.equal(unname(bright), unname(absolute))))
})

test_that('waveout is interpreted as log10(Angstrom)', {
  # At z = 0 the final wavelength grid is the rest-frame waveout grid.
  waveout = seq(2, 5.35, by = 0.01)
  sed = .ps_sed(z = 0, waveout = waveout)

  expect_equal(sed$FinalLum$wave, 10^waveout, tolerance = 1e-8)
  # Default grid spans 100 A to ~2.2e9 A, which is the mm/sub-mm regime.
  expect_equal(range(sed$FinalLum$wave), range(10^waveout), tolerance = 1e-8)
})

test_that('CHARLOT and FALL attenuation has the right pivot and normalisation', {
  # CF is exp(-tau * (wave/pivot)^pow), so at the pivot it is exp(-tau) and it
  # must be monotonically increasing with wavelength for pow < 0.
  expect_equal(CF(5500, tau = 0.3), exp(-0.3), tolerance = 1e-12)
  expect_equal(CF_birth(5500, tau = 1), exp(-1), tolerance = 1e-12)

  wave = seq(1000, 1e5, by = 500)
  expect_true(all(diff(CF(wave, tau = 1)) > 0))
  # More dust means more attenuation at every wavelength.
  expect_true(all(CF(wave, tau = 2) < CF(wave, tau = 1)))

  # The 2175 A bump can only add attenuation, never remove it.
  bump = CF_screen(c(2175.8, 3000), tau = 1, Eb = 1)
  no_bump = CF_screen(c(2175.8, 3000), tau = 1, Eb = 0)
  expect_true(all(bump < no_bump))
})

test_that('binlims returns contiguous bin edges spanning the input', {
  bins = binlims(c(1, 2, 4))

  expect_equal(bins$lo, c(0.5, 1.5, 3.0))
  expect_equal(bins$hi, c(1.5, 3.0, 5.0))
  # Edges meet, and the outer edges sit half a spacing beyond the data.
  expect_equal(bins$lo[-1], bins$hi[-length(bins$hi)])
})

test_that('greybody is a blackbody times a wavelength-dependent emissivity', {
  # 850 micron is 8.5e6 Angstrom, which is where the emissivity is unity.
  k850_wave = 8.5e6
  wave = 10^seq(3, 7.5, length.out = 800)

  bb = blackbody(wave, Temp = 50, k850 = 0.077)
  gb = greybody(wave, Temp = 50, beta = 1.5, k850 = 0.077)

  expect_equal(.k_lambda(k850_wave, 1.5), 1, tolerance = 1e-12)

  # k_lambda is (wave / 850um)^-beta, so it is *above* unity blueward of the
  # normalisation wavelength and below it to the red.
  expect_true(all(.k_lambda(wave[wave < k850_wave], 1.5) > 1))
  expect_true(all(.k_lambda(wave[wave > k850_wave], 1.5) < 1))

  # Because the emissivity cuts the blue side, the greybody peak sits blueward
  # of the blackbody peak (the Wien side is suppressed, the Rayleigh-Jeans side
  # is not).
  expect_lt(wave[which.max(gb)], wave[which.max(bb)])

  # A steeper beta suppresses the Wien side more, moving the peak further blue.
  peak = function(beta) wave[which.max(greybody(wave, Temp = 50, beta = beta))]
  expect_lt(peak(2.5), peak(1))
})

test_that('Inoue14_IGM is unity in the red and absorbing in the rest-frame UV', {
  # No absorption at z = 0.
  expect_equal(Inoue14_IGM(c(1215, 5000), 0), c(1, 1), tolerance = 1e-12)

  # At z = 3, rest-frame Lyman-alpha (1215 A) is strongly absorbed while a
  # redward wavelength is untouched.
  observed = c(1215, 5000) * 4
  trans = Inoue14_IGM(observed, 3)
  expect_lt(trans[1], 1)
  expect_equal(trans[2], 1, tolerance = 1e-10)
})