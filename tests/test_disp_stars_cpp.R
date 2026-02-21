# Regression test: verify that .disp_stars_cpp() matches the original R loop
# to within a tight numerical tolerance (relative error < 1e-10).
#
# Run via: Rscript tests/test_disp_stars_cpp.R
# (after the package has been installed or loaded with devtools::load_all())

library(ProSpect)

# R reference implementation of the disp_stars loop (identical to pre-Rcpp SFH.R)
disp_stars_R <- function(wave_lum_log, lum_log, z_disp, grid, weights, res) {
  lum_conv <- numeric(length(lum_log))
  for (i in seq_along(grid)) {
    z_seq_log  <- log10(1 + grid[i] * z_disp)
    new_wave   <- wave_lum_log + z_seq_log
    new_lum    <- lum_log - z_seq_log
    new_lum    <- 10^approx(x = new_wave, y = new_lum,
                            xout  = wave_lum_log,
                            rule  = 2,
                            yleft = new_lum[1],
                            yright = new_lum[length(new_lum)])$y
    lum_conv <- lum_conv + weights[i] * new_lum
  }
  lum_conv * res
}

set.seed(42)
n    <- 500
wave <- seq(3000, 10000, length.out = n)
wave_log <- log10(wave)
lum      <- runif(n, 1e3, 1e8)
lum_log  <- log10(lum)

# Several test cases with different dispersion and grid settings
test_cases <- list(
  list(veldisp = 100,  range = 3, res = 0.1),
  list(veldisp = 250,  range = 3, res = 0.1),
  list(veldisp = 50,   range = 5, res = 0.05),
  list(veldisp = 300,  range = 3, res = 0.2)
)

c_to_mps <- 299792458  # speed of light in m/s (matches ProSpect:::.c_to_mps)

all_passed <- TRUE

for (tc in test_cases) {
  veldisp <- tc$veldisp
  range   <- tc$range
  res     <- tc$res

  grid    <- seq(-range, range, by = res)
  weights <- dnorm(grid)
  z_disp  <- rep(veldisp / (c_to_mps / 1000), n)  # constant z_disp for simplicity

  ref <- disp_stars_R(wave_log, lum_log, z_disp, grid, weights, res)
  got <- .disp_stars_cpp(wave_log, lum_log, z_disp, grid, weights, res)

  rel_err <- max(abs(ref - got) / pmax(abs(ref), 1e-300))

  if (rel_err > 1e-10) {
    cat(sprintf("FAIL: veldisp=%g range=%g res=%g => max rel error = %e\n",
                veldisp, range, res, rel_err))
    all_passed <- FALSE
  } else {
    cat(sprintf("PASS: veldisp=%g range=%g res=%g => max rel error = %e\n",
                veldisp, range, res, rel_err))
  }
}

if (!all_passed) stop("One or more disp_stars_cpp regression tests failed!")
cat("All disp_stars_cpp regression tests passed.\n")
