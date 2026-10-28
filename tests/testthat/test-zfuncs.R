# Tests for the mass-mapping metallicity histories, Zfunc_massmap_lin and
# Zfunc_massmap_box.
#
# Both build a cumulative mass distribution on a fixed 1e8-year grid. The
# normalisation used to come from adaptive stats::integrate() while the CDF
# itself was a rectangle sum on that grid, so the two disagreed by ~0.8% and the
# CDF was clipped at 0 to hide it. integrate() could also miss a narrow burst
# entirely and return 0 with message "OK", which silently collapsed the whole
# history to a constant Zfinal. These tests pin the corrected behaviour.

.ps_zages = function() {
  .ps_speclib()$Age
}

# A burst far narrower than the 1e8-year CDF grid, sitting on a grid point so
# it cannot be missed by sampling: finite, non-negative, and huge in amplitude.
.ps_narrow_burst = function(centre = 5e9, width = 1e6, amp = 1e6) {
  force(centre); force(width); force(amp)
  function(age, ...) amp * exp(-((age - centre) / width)^2)
}

test_that('mass-mapped histories vary with age rather than collapsing', {
  ages = .ps_zages()

  # A smooth SFH must give a metallicity that actually evolves.
  for (massfunc in list(massfunc_b5, massfunc_snorm, massfunc_exp)) {
    Z = Zfunc_massmap_lin(ages, massfunc = massfunc)

    expect_false(anyNA(Z))
    expect_gt(diff(range(Z)), 1e-3)
    expect_true(all(Z >= 1e-4))
    expect_true(all(Z <= 0.02 + 1e-9))
  }
})

test_that('a narrow burst is not silently zeroed away', {
  ages = .ps_zages()

  # This is the regression: integrate() returned 0 here, so masssum > 0 was
  # FALSE and the function fell through to rep(Zfinal, ...).
  Z = Zfunc_massmap_lin(ages, massfunc = .ps_narrow_burst())

  expect_false(anyNA(Z))
  # A burst is a step in cumulative mass, so the history must be bimodal
  # between the floor and Zfinal rather than flat.
  expect_gt(length(unique(round(Z, 8))), 1)
  expect_gt(diff(range(Z)), 1e-3)

  # The same holds for the box mapping.
  Zbox = Zfunc_massmap_box(ages, massfunc = .ps_narrow_burst())
  expect_false(anyNA(Zbox))
  expect_gt(diff(range(Zbox)), 1e-3)
})

test_that('Zfunc_massmap_lin endpoints are Zfinal today and Zstart at the start', {
  ages = .ps_zages()

  Z = Zfunc_massmap_lin(ages, massfunc = massfunc_b5)

  # Youngest population has the most metals, oldest the least.
  expect_gt(Z[which.min(ages)], Z[which.max(ages)])
  expect_equal(Z[which.max(ages)], 1e-4, tolerance = 1e-6)
  # And the youngest bin is close to Zfinal for a b5 history.
  expect_gt(max(Z), 0.015)
})

test_that('ages beyond Zagemax are pinned to Zstart', {
  ages = .ps_zages()

  Z = Zfunc_massmap_lin(ages, massfunc = massfunc_b5, Zagemax = 5)
  expect_true(all(Z[ages > 5e9] == 1e-4))

  Zbox = Zfunc_massmap_box(ages, massfunc = massfunc_b5, Zagemax = 5)
  expect_true(all(Zbox[ages > 5e9] == 1e-4))
})

test_that('Zagemax is inferred from magemax when not given', {
  ages = .ps_zages()

  # massfunc_* take magemax in Gyr; Zfunc_massmap_* should honour it via ...
  Z = Zfunc_massmap_lin(ages, massfunc = massfunc_b5, magemax = 5)
  expect_true(all(Z[ages > 5e9] == 1e-4))
})

test_that('Zstart and Zfinal bracket the returned history', {
  ages = .ps_zages()

  Z = Zfunc_massmap_lin(ages, massfunc = massfunc_b5, Zstart = 0.005, Zfinal = 0.03)
  expect_gte(min(Z), 0.005 - 1e-9)
  expect_lte(max(Z), 0.03 + 1e-9)

  Zbox = Zfunc_massmap_box(ages, massfunc = massfunc_b5, Zstart = 0.005, Zfinal = 0.03)
  expect_gte(min(Zbox), 0.005 - 1e-9)
})

test_that('a zero SFH falls back to a constant Zfinal', {
  ages = .ps_zages()
  zero = function(age, ...) rep(0, length(age))

  expect_equal(Zfunc_massmap_lin(ages, massfunc = zero), rep(0.02, length(ages)))
  expect_equal(Zfunc_massmap_box(ages, massfunc = zero), rep(0.02, length(ages)))
})

test_that('massfunc is required', {
  expect_error(Zfunc_massmap_lin(.ps_zages()), 'Need massfunc')
  expect_error(Zfunc_massmap_box(.ps_zages()), 'Need massfunc')
})

test_that('a higher Zfinal gives a higher metallicity everywhere', {
  ages = .ps_zages()

  low = Zfunc_massmap_lin(ages, massfunc = massfunc_b5, Zfinal = 0.01)
  high = Zfunc_massmap_lin(ages, massfunc = massfunc_b5, Zfinal = 0.03)

  expect_true(all(high >= low))
  expect_gt(mean(high), mean(low))
})

test_that('the box mapping predicts a lower Z than the linear one at fixed Zfinal', {
  ages = .ps_zages()

  lin = Zfunc_massmap_lin(ages, massfunc = massfunc_b5)
  box = Zfunc_massmap_box(ages, massfunc = massfunc_b5)

  # The box mapping is logarithmic in gas fraction, so it is bounded by
  # Zfinal as the gas runs out and should sit below the linear map mid-history.
  expect_true(all(is.finite(box)))
  expect_lt(mean(box), mean(lin))
})

test_that('the Z functions work inside ProSpectSED', {
  # End-to-end: the history must actually reach the SED model as a varying Z.
  sed = .ps_sed(Z = Zfunc_massmap_lin)

  expect_true(is.finite(sed$Stars$masstot))
  expect_gt(sed$Stars$masstot, 0)
  expect_gt(diff(range(sed$Stars$Zvec)), 1e-3)
  expect_true(all(is.finite(sed$Photom)))
})