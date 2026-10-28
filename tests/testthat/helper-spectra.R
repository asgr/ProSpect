# Shared fixtures for the ProSpect test suite.
#
# Everything here is deliberately cheap: `sparse = 10` on the BC03lr library
# gives a model that runs in well under a second and is accurate enough for
# qualitative/physical assertions.

# 17 bands spanning the FUV to the sub-mm. The FIR bands matter because they are
# where dust re-emission shows up.
.ps_bands = c(
  'FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'z_SDSS',
  'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE',
  'P100_Herschel', 'S250_Herschel', 'S500_Herschel'
)

# Bands dominated by starlight, which dust can only suppress.
.ps_stellar_bands = .ps_bands[1:11]

# Bands dominated by dust re-emission, which brighten as optical depth rises.
.ps_dust_bands = c('W3_WISE', 'W4_WISE', 'P100_Herschel', 'S250_Herschel', 'S500_Herschel')

.ps_filtout = function(bands = .ps_bands) {
  out = lapply(bands, function(band) approxfun(getfilt(band)))
  names(out) = bands
  out
}

# Default ProSpectSED arguments, overridden by anything passed through ...
# Uses whole-value replacement rather than modifyList(), which would recurse
# into list arguments like filtout instead of substituting them.
.ps_args = function(...) {
  user = list(...)
  if (length(user) > 0 && (is.null(names(user)) || any(names(user) == ''))) {
    stop('All arguments passed to the ProSpect test helpers must be named.')
  }

  args = list(
    massfunc = massfunc_b5,
    speclib = .ps_speclib(),
    Dale = .ps_dale(),
    filtout = .ps_filtout(),
    filters = NULL,
    z = 0.1,
    sparse = 10
  )

  for (name in names(user)) {
    args[[name]] = user[[name]]
  }

  args
}

# Generate an SED and return Photom as a vector named by filter.
.ps_photom = function(...) {
  args = .ps_args(...)
  out = do.call(ProSpectSED, args)
  stats::setNames(out$Photom, names(args$filtout))
}

# Full ProSpectSED list, for tests that need Stars / DustEmit / dustlum etc.
.ps_sed = function(...) {
  do.call(ProSpectSED, .ps_args(...))
}

.ps_speclib = function() {
  if (is.null(.ps_cache$speclib)) {
    e = new.env()
    data('BC03lr', envir = e)
    .ps_cache$speclib = e$BC03lr
  }
  .ps_cache$speclib
}

.ps_dale = function() {
  if (is.null(.ps_cache$dale)) {
    e = new.env()
    data('Dale_NormTot', envir = e)
    .ps_cache$dale = e$Dale_NormTot
  }
  .ps_cache$dale
}

.ps_cache = new.env(parent = emptyenv())

# Integrate a wave/lum data frame with the same left-rectangle rule the package
# uses internally (.qdiff), so energy-budget tests compare like with like.
.ps_integral = function(spec) {
  wave = spec[[1]]
  lum = spec[[2]]
  sum(c(0, diff(wave)) * lum)
}

# Ratio of two bands in a named photometry vector, as a plain numeric.
.ps_ratio = function(photom, num, den) {
  unname(photom[num] / photom[den])
}