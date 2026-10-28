# Argument-validation tests for filterTranMags(). The mag_out check used to be
# `if(!is.vector(mag_out)){'mag_out should be a vector!'}` -- a bare string,
# which evaluated and discarded, so bad input was silently accepted.

.ps_mock_mags = function(n = 200) {
  set.seed(9114)
  g = stats::runif(n, 18, 24)
  # r correlates with g plus a colour term, so the fit is well determined.
  r = g - 0.5 - 0.3 * stats::rnorm(n, 0, 0.4)
  data.frame(g_SDSS = g, r_SDSS = r)
}

# A third band with its own scatter. It must not be a deterministic function of
# r, or the colour term is constant and lm() cannot fit it (returning NA).
.ps_mock_i = function(mags) {
  set.seed(4420)
  mags$r_SDSS - 0.1 + stats::rnorm(nrow(mags), 0, 0.15)
}

test_that('filterTranMags rejects a non-vector mag_out', {
  mags = .ps_mock_mags()

  # Every one of these used to pass through silently.
  expect_error(filterTranMags(mags$g_SDSS, mags), 'mag_out must be a vector')
  expect_error(filterTranMags(mags$g_SDSS, as.matrix(mags)), 'mag_out must be a vector')
})

test_that('filterTranMags still accepts a plain vector mag_out', {
  mags = .ps_mock_mags()

  single = filterTranMags(mags$g_SDSS, mags$r_SDSS)
  expect_type(single, 'double')
  expect_named(single, c('beta', 'sigma'))

  # The recovered offset should be close to the -0.5 used to build the mock.
  expect_equal(unname(single['beta']), -0.5, tolerance = 0.1)
})

test_that('filterTranMags recovers a colour term from multi-band input', {
  mags = .ps_mock_mags()
  mags$i_SDSS = .ps_mock_i(mags)

  out = filterTranMags(as.matrix(mags), mags$r_SDSS)
  expect_type(out, 'list')
  expect_gt(length(out), 0)

  # Each entry is a named alpha/beta/sigma triple.
  for (fit in out) {
    expect_named(fit, c('alpha', 'beta', 'sigma'))
    expect_true(all(is.finite(fit)))
  }
})

test_that('filterTranMags rejects an unknown return option', {
  mags = .ps_mock_mags()
  mags$i_SDSS = mags$r_SDSS

  expect_error(
    filterTranMags(as.matrix(mags), mags$r_SDSS, return = 'not-an-option'),
    'return must be one of'
  )
})

test_that('filterTranMags best-fit selectors return the prediction contract', {
  mags = .ps_mock_mags()
  mags$i_SDSS = .ps_mock_i(mags)

  for (option in c('bestalpha', 'bestbeta', 'bestsigma', 'bestall')) {
    out = filterTranMags(as.matrix(mags), mags$r_SDSS, return = option)

    expect_named(out, c('params', 'ref_bands', 'predict'),
                 label = paste(option, 'returns the documented list'))
    expect_length(out$predict, nrow(mags))
    expect_true(all(is.finite(out$predict)), label = paste(option, 'predictions finite'))
    expect_length(out$ref_bands, 3)
  }
})
# --- filterTranBands -------------------------------------------------------
#
# filterTranBands samples random SFHs, converts them to magnitudes and hands
# them to filterTranMags. It builds a column name for the target band with
# names('Target') -- which is always NULL -- so an unnamed filt_out produced an
# unnamed (actually NA-named) target column. Nsamp is tiny here: the point is
# the plumbing, not the statistics.

.ps_two_bands = function() {
  list(g_SDSS = approxfun(getfilt('g_SDSS')),
       r_SDSS = approxfun(getfilt('r_SDSS')))
}

test_that('the target column name falls back to Target when filt_out is unnamed', {
  # A bare function is the documented and usual way to pass filt_out.
  expect_equal(ProSpect:::.tranTargetName(approxfun(getfilt('r_VST'))), 'Target')

  # names('Target') is NULL, which is what caused the bug.
  expect_null(names('Target'))
})

test_that('the target column name is used when filt_out is a named list', {
  expect_equal(ProSpect:::.tranTargetName(list(r_VST = approxfun(getfilt('r_VST')))),
               'r_VST')

  # Anything other than exactly one name falls back rather than partially naming.
  expect_equal(ProSpect:::.tranTargetName(list(approxfun(getfilt('r_VST')))), 'Target')
  expect_equal(ProSpect:::.tranTargetName(list(a = 1, b = 2)), 'Target')
})

test_that('filterTranBands runs for the documented input shapes', {
  # Multiple named bands against a bare target function. This is the shape used
  # in the man page example, and the one that previously produced an unnamed
  # target column.
  multi = filterTranBands(.ps_two_bands(), approxfun(getfilt('r_VST')), Nsamp = 3)

  expect_type(multi, 'list')
  expect_gt(length(multi), 0)
  # Every returned mapping names real input bands, not a placeholder NA.
  expect_true(all(grepl('g_SDSS|r_SDSS', names(multi))))

  # A single bare input function reduces to an offset-only fit.
  single = filterTranBands(approxfun(getfilt('g_SDSS')), approxfun(getfilt('r_VST')),
                           Nsamp = 3)
  expect_named(single, c('beta', 'sigma'))
  expect_true(all(is.finite(single)))
})

test_that('filterTranBands still demands names for a multi-band filt_in', {
  expect_error(
    filterTranBands(list(approxfun(getfilt('g_SDSS')), approxfun(getfilt('r_SDSS'))),
                    approxfun(getfilt('r_VST')), Nsamp = 3),
    'filt_in must have filter names'
  )
})
