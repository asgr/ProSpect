# Tests for ProSpectSEDlike(), the likelihood interface used by optim, CMA-ES,
# LaplacesDemon and Highlander. The key property is consistency: the likelihood
# must reproduce exactly what ProSpectSED generates at the same parameters, and
# must be maximal there. If those two drift apart, every fit silently
# optimises to the wrong answer.

.ps_like_bands = c('g_SDSS', 'r_SDSS', 'K_VISTA', 'W3_WISE')

# Build a mock catalogue at known parameters plus a matching Data list.
# Everything the likelihood needs must be pinned explicitly, because it falls
# back to defaults different from this package's generating defaults.
.ps_mock_fit = function(truth = c(tau_birth = 0.5, tau_screen = 0.3),
                        noise = 0.05) {
  filtout = .ps_filtout(.ps_like_bands)
  LumDist_Mpc = celestial::cosdistLumDist(z = 0.1, H0 = 67.8, OmegaM = 0.308)

  args = .ps_args(
    filtout = filtout,
    sparse = 5,
    LumDist_Mpc = LumDist_Mpc,
    tau_birth = unname(truth['tau_birth']),
    tau_screen = unname(truth['tau_screen'])
  )
  flux = do.call(ProSpectSED, args)$Photom

  list(
    truth = truth,
    flux = data.frame(filter = .ps_like_bands, pivwave = 1, flux = flux,
                      fluxerr = flux * noise),
    Data = list(
      flux = data.frame(filter = .ps_like_bands, pivwave = 1, flux = flux,
                        fluxerr = flux * noise),
      filtout = filtout,
      parm.names = names(truth),
      speclib = .ps_speclib(),
      Dale = .ps_dale(),
      SFH = SFHfunc,
      arglist = list(z = 0.1, massfunc = massfunc_b5, LumDist_Mpc = LumDist_Mpc)
    )
  )
}

# Call the likelihood with a given fit mode and parameter vector.
.ps_like = function(mock, parm = mock$truth, fit = 'optim', ...) {
  ProSpectSEDlike(parm, utils::modifyList(mock$Data, c(list(fit = fit), list(...))))
}

test_that('the likelihood reproduces the generating SED exactly at the truth', {
  mock = .ps_mock_fit()

  # fit = 'check' returns the full model, so we can compare photometry directly.
  checked = .ps_like(mock, fit = 'check')

  expect_s3_class(checked, 'ProSpectSEDlike')
  expect_equal(unname(checked$SEDout$Photom), mock$flux$flux, tolerance = 1e-12)
  # ...and therefore residuals of exactly zero.
  expect_lt(max(abs(checked$SEDout$Photom / mock$flux$flux - 1)), 1e-12)

  # Dev is -2*log-likelihood (not chi-squared), so it is only finite, not zero.
  expect_true(is.finite(checked$Dev))
  expect_true(is.finite(checked$LP))
})

test_that('the likelihood is maximised at the generating parameters', {
  mock = .ps_mock_fit()

  # optim/cma return -LP, so a worse fit means a larger value.
  at_truth = .ps_like(mock)
  perturbed = c(
    .ps_like(mock, parm = c(tau_birth = 0.3, tau_screen = 0.3)),
    .ps_like(mock, parm = c(tau_birth = 1.0, tau_screen = 0.3)),
    .ps_like(mock, parm = c(tau_birth = 0.5, tau_screen = 1.0)),
    .ps_like(mock, parm = c(tau_birth = 0.0, tau_screen = 0.0))
  )

  expect_true(all(perturbed > at_truth))
  # As a minimum, the objective must be finite and scalar for an optimiser.
  expect_length(at_truth, 1)
  expect_true(is.finite(at_truth))
})

test_that('the likelihood returns the shape each optimiser expects', {
  mock = .ps_mock_fit()

  # optim / CMA-ES want a scalar to minimise
  expect_type(.ps_like(mock, fit = 'optim'), 'double')
  expect_type(.ps_like(mock, fit = 'cma'), 'double')

  # LaplacesDemon / Highlander want LP, Dev, Monitor, yhat, parm
  ld = .ps_like(mock, fit = 'ld')
  expect_type(ld, 'list')
  expect_named(ld, c('LP', 'Dev', 'Monitor', 'yhat', 'parm'))
  expect_equal(ld$Dev, -2 * ld$LP, tolerance = 1e-8)
  expect_equal(unname(ld$parm), unname(mock$truth))
})

test_that('Monitor only carries the quantities requested in mon.names', {
  mock = .ps_mock_fit()

  # masstot is a derived quantity, so it must be reported explicitly...
  monitored = .ps_like(mock, fit = 'ld', mon.names = 'masstot')
  expect_true('masstot' %in% names(monitored$Monitor))
  expect_gt(monitored$Monitor['masstot'], 0)

  # ...and an empty request still returns something non-NULL, since LD cannot
  # handle a zero-length Monitor.
  unmonitored = .ps_like(mock, fit = 'ld', mon.names = character(0))
  expect_false(is.null(unmonitored$Monitor))
})

test_that('the likelihood rejects non-scalar fitted parameters', {
  mock = .ps_mock_fit()

  # Non-scalar values in arglist reach ProSpectSED as vectorised arguments and
  # fail somewhere inside the model rather than at a clean guard. This test
  # documents the current behaviour; a proper stop('... must be scalar')
  # would be an improvement.
  expect_error(
    suppressWarnings(
      .ps_like(mock, arglist = list(z = c(0.1, 0.2), massfunc = massfunc_b5,
                                    LumDist_Mpc = 400))
    )
  )
})

test_that('likelihood priors and constraints are applied', {
  mock = .ps_mock_fit()

  no_prior = .ps_like(mock, fit = 'ld')
  with_prior = .ps_like(mock, fit = 'ld',
                        prior = function(parm) -sum((parm - 0.5)^2))
  expect_lt(with_prior$LP, no_prior$LP)

  # A constraint mapping is applied to the parameters before the model is built.
  with_constraint = .ps_like(mock, fit = 'ld',
                             constraints = function(parm) { parm['tau_screen'] = 1; parm })
  expect_gt(with_constraint$Dev, 1e-8)
})

test_that('a bad like option is a hard error', {
  mock = .ps_mock_fit()
  expect_error(.ps_like(mock, like = 'not-a-likelihood'), 'Bad like option')
})

test_that('the three likelihood families all give a finite objective', {
  mock = .ps_mock_fit()

  for (like in c('st', 'norm', 'chisq')) {
    value = .ps_like(mock, like = like)
    expect_true(is.finite(value), label = paste('finite objective for', like))
  }
})