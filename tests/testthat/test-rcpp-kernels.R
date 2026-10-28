# Regression tests for the Rcpp kernels against plain-R reference
# implementations. These guard the numerical equivalence of the C++ fast paths:
# if a kernel is recompiled with a bug, the R reference stays correct and the
# comparison fails.

# Reference implementation of .disp_stars_cpp: the loop that used to live in
# SFHfunc()/SFHburst() before the kernel was introduced.
disp_stars_R = function(wave_lum_log, lum_log, z_disp, grid, weights, res) {
  lum_conv = numeric(length(lum_log))
  for (i in seq_along(grid)) {
    z_seq_log = log10(1 + grid[i] * z_disp)
    new_wave = wave_lum_log + z_seq_log
    new_lum = lum_log - z_seq_log
    new_lum = 10^stats::approx(x = new_wave, y = new_lum,
                               xout = wave_lum_log,
                               rule = 2,
                               yleft = new_lum[1],
                               yright = new_lum[length(new_lum)])$y
    lum_conv = lum_conv + weights[i] * new_lum
  }
  lum_conv * res
}

test_that('.disp_stars_cpp matches the R reference across grid settings', {
  set.seed(3891)
  n = 500
  wave_log = log10(seq(3000, 10000, length.out = n))
  lum_log = log10(stats::runif(n, 1e3, 1e8))

  cases = list(
    list(veldisp = 100, range = 3, res = 0.1),
    list(veldisp = 250, range = 3, res = 0.1),
    list(veldisp = 50, range = 5, res = 0.05),
    list(veldisp = 300, range = 3, res = 0.2)
  )

  for (case in cases) {
    grid = seq(-case$range, case$range, by = case$res)
    weights = stats::dnorm(grid)
    z_disp = rep(case$veldisp / (299792458 / 1000), n)

    ref = disp_stars_R(wave_log, lum_log, z_disp, grid, weights, case$res)
    got = ProSpect:::.disp_stars_cpp(wave_log, lum_log, z_disp, grid, weights, case$res)

    rel_err = max(abs(ref - got) / pmax(abs(ref), 1e-300))
    expect_lt(rel_err, 1e-10)
  }
})

test_that('.disp_stars_cpp conserves total flux for a zero grid', {
  # With only the zero-velocity grid point the convolution is the identity.
  wave_log = log10(seq(3000, 10000, length.out = 200))
  lum_log = log10(seq(1e3, 1e6, length.out = 200))

  got = ProSpect:::.disp_stars_cpp(wave_log, lum_log, rep(0.1, 200), 0, 1, 1)
  expect_equal(got, 10^lum_log, tolerance = 1e-12)
})

test_that('.spec_rebin_cpp conserves flux when rebinning a flat spectrum', {
  # A constant f_lambda must stay constant under rebinning, to within the
  # discretisation of the input grid. The kernel returns a bare numeric vector
  # (no names), so index it directly.
  wave = seq(3000, 10000, by = 1)
  flux = rep(1.5, length(wave))
  grid = seq(3500, 9000, by = 250)

  out = ProSpect:::.spec_rebin_cpp(wave, flux, grid, logbin_in = FALSE, logbin_out = FALSE)
  expect_length(out, length(grid))
  expect_equal(as.numeric(out), rep(1.5, length(grid)), tolerance = 1e-8)
})

test_that('.spec_rebin_cpp returns flux and weight columns when invar is given', {
  wave = seq(3000, 10000, by = 1)
  flux = rep(1.5, length(wave))
  grid = seq(3500, 9000, by = 250)
  invar = rep(1, length(wave))

  out = ProSpect:::.spec_rebin_cpp(wave, flux, grid, invar_in = invar,
                                   logbin_in = FALSE, logbin_out = FALSE)
  # flux in the first half, accumulated inverse variance in the second
  expect_length(out, 2 * length(grid))
  expect_equal(as.numeric(out[seq_along(grid)]), rep(1.5, length(grid)), tolerance = 1e-8)
  expect_true(all(out[length(grid) + seq_along(grid)] > 0))
})

test_that('.vec_add_cpp and .vec_mult_cpp match base vectorised arithmetic', {
  # Both kernels mutate their first argument in place and return invisibly.
  # They do NOT copy, so the caller's object is modified too -- package code
  # only ever passes freshly-built vectors. Always pass a fresh vector here.
  set.seed(7723)
  a = stats::runif(1000, 1, 100)
  b = stats::runif(1000, 1, 100)

  added = a + 0   # fresh copy, because the kernel mutates in place
  ProSpect:::.vec_add_cpp(added, b)
  expect_equal(added, a + b, tolerance = 1e-14)

  scaled = a + 0
  ProSpect:::.vec_mult_cpp(scaled, b)
  expect_equal(scaled, a * b, tolerance = 1e-14)

  # Scalar second argument is broadcast.
  plus_two = a + 0
  ProSpect:::.vec_add_cpp(plus_two, 2)
  expect_equal(plus_two, a + 2, tolerance = 1e-14)

  # Mismatched lengths are a hard error.
  expect_error(ProSpect:::.vec_add_cpp(a + 0, b[1:3]))
  expect_error(ProSpect:::.vec_mult_cpp(a + 0, b[1:3]))
})

test_that('the in-place kernels mutate their input (documented hazard)', {
  # Locking this down so an unexpected copy is noticed: if either kernel ever
  # gains copy semantics, package code relying on in-place updates changes
  # behaviour silently.
  x = c(1, 2, 3)
  ProSpect:::.vec_add_cpp(x, c(10, 20, 30))
  expect_equal(x, c(11, 22, 33))

  y = c(1, 2, 3)
  ProSpect:::.vec_mult_cpp(y, 2)
  expect_equal(y, c(2, 4, 6))
})

test_that('.colSums_wt_cpp matches weighted colSums', {
  set.seed(4410)
  mat = matrix(stats::runif(200 * 30, 1, 10), nrow = 200)
  wt = stats::runif(200, 0.1, 2)

  expect_equal(as.numeric(ProSpect:::.colSums_wt_cpp(mat, wt)),
               as.numeric(crossprod(mat, wt)),
               tolerance = 1e-10)
  # Unit weight is a plain colSums.
  expect_equal(as.numeric(ProSpect:::.colSums_wt_cpp(mat)),
               as.numeric(colSums(mat)),
               tolerance = 1e-10)
})