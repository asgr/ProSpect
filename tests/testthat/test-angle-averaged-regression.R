# Standalone: Rscript --vanilla -e \
#   'testthat::test_file("tests/testthat/test-angle-averaged-regression.R")'
# For namespace tests, explicitly load the target with pkgload::load_all(path)
# and pass package = 'ProSpect' to test_file(). R CMD check selects its installed
# target automatically.
.aa_source = function(root) {
  env = new.env(parent = asNamespace('stats'))
  env$data = function(..., package = 'ProSpectData') {
    utils::data(..., package = package)
  }
  for (name in c('constants', 'utility', 'photom', 'dust', 'AGN',
                 'radiocont', 'ProSpect')) {
    sys.source(file.path(root, 'R', paste0(name, '.R')), envir = env)
  }
  env
}
.aa_code = local({
  # Only trust the namespace selected by the package test runner. A standalone
  # test_file() must not silently test a pre-existing global installation.
  if (identical(testthat::testing_package(), 'ProSpect')) {
    if (!'ProSpect' %in% loadedNamespaces()) {
      stop('Load the target ProSpect with pkgload::load_all(path) before testing.')
    }
    return(asNamespace('ProSpect'))
  }

  # Standalone source-tree runs remain isolated from any loaded ProSpect.
  path = normalizePath(getwd(), mustWork = TRUE)
  repeat {
    description = file.path(path, 'DESCRIPTION')
    if (file.exists(file.path(path, 'R', 'ProSpect.R')) &&
        file.exists(description) &&
        identical(unname(read.dcf(description, fields = 'Package')[1, 1]),
                  'ProSpect')) break
    parent = dirname(path)
    if (identical(parent, path)) {
      stop(paste('Run inside the ProSpect source tree, or use',
                 'pkgload::load_all(target_path) and',
                 "testthat::test_file(test_path, package = 'ProSpect')."))
    }
    path = parent
  }
  .aa_source(path)
})

# Tiny deterministic grids use the upstream interpolators' native dimensions.
# Mean and viewing-angle wavelengths deliberately differ.
.aa_templates = function(model = 'Fritz') {
  wave = 10^seq(2, 10, length.out = if (model == 'Fritz') 178 else 132)
  meanwave = 10^seq(2.1, 9.9, length.out = 93)
  view_shape = (wave / 5500)^-1.3 * exp(-wave / 2e6) +
    0.003 * (wave / 1e6)^-0.5
  mean_shape = 2.7 * ((meanwave / 5500)^-0.9 * exp(-meanwave / 8e5) +
    0.004 * (meanwave / 1e6)^-0.5)
  if (model == 'Fritz') {
    axes = list(ct = c(20, 60), al = c(0, 6), be = c(-1, 0),
                ta = c(0.1, 6), rm = c(10, 100))
    view = c(axes, list(an = c(0, 90), Wave = wave,
      Aspec = array(rep(view_shape * 1e-5, each = 64) * rep(c(1, 3), each = 32),
                    c(rep(2, 6), 178))))
    mean = c(axes, list(Wave = meanwave,
      Aspec = array(rep(mean_shape * 1e-5, each = 32), c(rep(2, 5), 93))))
  } else {
    axes = list(ta = c(0.1, 6), p = c(0, 2), q = c(0, 2),
                ct = c(20, 60), rm = c(10, 100))
    view = c(axes, list(an = c(0, 90), Wave = wave,
      Aspec = array(rep(view_shape * wave * 3.839e39, each = 64) *
                      rep(c(1, 3), each = 32),
                    c(rep(2, 5), 1, 2, 132))))
    mean = c(axes, list(Wave = meanwave,
      Aspec = array(rep(mean_shape * meanwave * 3.839e39, each = 32),
                    c(rep(2, 5), 1, 93))))
  }
  class(view) = model
  class(mean) = paste0(model, '_anav')
  list(view = view, mean = mean)
}

.aa_dale = function() {
  wave = 10^seq(2, 10, length.out = 1496)
  shape = exp(-0.5 * ((log10(wave) - 6) / 0.4)^2) / wave + 1e-18
  shape = shape / sum(shape * c(0, diff(wave)))
  list(Wave = wave, alpha_SF = c(0.0625, 4), AGNfrac = 0,
       Aspec = list(rbind(shape, shape)))
}

.aa_stars = function(tau_birth, tau_screen, pow_birth, pow_screen, ...) {
  wave = 10^seq(2, 10, length.out = 211)
  unatten = 1e5 * (wave / 5500)^-1.5 * exp(-wave / 1e5) + 1e-12
  birth = unatten * (0.4 + 0.6 * exp(-tau_birth * (wave / 5500)^pow_birth))
  atten = birth * exp(-tau_screen * (wave / 5500)^pow_screen)
  integral = function(x) sum(x * c(0, diff(wave)))
  list(wave = wave, wave_lum = wave, lum_unatten = unatten,
       lum_atten = atten, lumtot_birth = integral(unatten - birth),
       lumtot_screen = integral(birth - atten), masstot = 1e10)
}

.aa_args = function(model = 'Fritz') {
  template = .aa_templates(model)
  list(SFH = .aa_stars, AGN = template$view, AGNlum = 1e44,
       Dale = .aa_dale(), Dale_M2L_func = function(alpha_SF) 2,
       z = 0.1, LumDist_Mpc = 100, filters = NULL,
       waveout = c(2, 3, 4, 5, 6, 7, 8, 9, 9.35), returnall = TRUE)
}

.aa_none_args = function(model, grid) {
  args = .aa_args(if (model == 'legacy') 'Fritz' else model)
  if (model == 'legacy') {
    args$AGN = list(Wave = 10^seq(2, 10, length.out = 178),
                    Aspec = rep(1e-9, 178))
  }
  args$tau_birth = 0.6
  args$tau_screen = 0.8
  args$pow_screen = -1.1
  args$tau_AGN = 1.4
  args$pow_AGN = -0.9
  args$alpha_SF_AGN = 0.7
  args$alpha_SF_screen = 2
  args$AGNan = 65
  args$Eb = 1.5
  args$L0 = 2200
  args$LFWHM = 600
  args$ff_frac_SF = 0.18
  args$ff_frac_AGN = 0.23
  args$Te_AGN = 15000
  args$sy_power_AGN = -0.95
  # Neither None nor viewing_angle may inspect or load a mean template.
  args$AGN_anav = list(unused = TRUE)
  if (grid == 'NULL') args['waveout'] = list(NULL)
  args
}

.aa_run = function(args, ..., code = .aa_code) {
  overrides = list(...)
  args[names(overrides)] = overrides
  do.call(code$ProSpectSED, args)
}

# Rebind copies of the entry point and helper so the data spy never modifies a
# loaded package namespace (and works for standalone source-tree tests too).
.aa_loader_fixture = function(model, loaded = .aa_templates(model)$mean) {
  code = new.env(parent = .aa_code)
  calls = new.env(parent = emptyenv())
  calls$names = character()
  calls$from_helper = logical()
  helper_name = paste0(model, '_anav_interp')
  code$data = function(name, envir, ...) {
    calls$names = c(calls$names, name)
    calls$from_helper = c(calls$from_helper,
                          identical(sys.function(-1), code[[helper_name]]))
    assign(name, loaded, envir = envir)
    invisible(name)
  }
  for (name in c('ProSpectSED', helper_name)) {
    fun = .aa_code[[name]]
    environment(fun) = code
    code[[name]] = fun
  }
  list(code = code, calls = calls)
}

.aa_error = function(expr, pattern) {
  error = tryCatch(force(expr), error = identity)
  testthat::expect_s3_class(error, 'error')
  if (inherits(error, 'error')) {
    testthat::expect_match(conditionMessage(error), pattern, ignore.case = TRUE)
  }
}

.aa_on_grid = function(wave, flux, output_wave) {
  10^stats::approx(log10(wave), log10(flux), xout = log10(output_wave),
                  yleft = -Inf, yright = -Inf)$y
}

.aa_finite = function(out) {
  testthat::expect_s3_class(out, 'ProSpectSED')
  for (name in c('FinalLum', 'AGN', 'Photom')) {
    spec = out[[name]]
    if (is.null(spec)) next
    testthat::expect_identical(all(is.finite(as.matrix(spec))), TRUE,
                               info = name)
    testthat::expect_gt(nrow(spec), 0L)
  }
}

# Golden values generated from commit dcfac56cc38873be697e571d00f8b7faf58440db,
# sourcing its R files with .aa_source() and running .aa_args(). These constants
# survive deletion of the baseline checkout and need no installed ProSpect/data.
.aa_golden = list(
  Fritz = list(
    final = c(674803.322342779, 1699351.70436128, 195710.911159956,
              17391.1081207665, 10299.960672994, 11.8191908844116,
              0.000100251404459963, 2.69273424634253e-9, 0.92784765882147),
    agn = c(560358.140569772, 1498950.08272078, 174184.189255069,
            15786.9968825434, 7609.92717977816, 2.99792458e-196,
            2.99792458e-198, 2.99792458e-200, 0.927847656130377),
    dustlum = c(3949353348.49075, 16016046789.5488, NA, 19965400138.0396),
    photom = c(7.92153760490617e-7, 0.000199487435588251, 0.00229745659381404,
               0.0204154769854054, 1.20911565042252, 0.138745856682181,
               0.000117685441680569, 3.16100927280049e-7, 545.89479365032)
  ),
  SKIRTOR = list(
    final = c(674803.32234278, 1698120.07925775, 195700.54151026,
              17318.512472251, 10135.5277799897, 11.8191908844116,
              0.000100251404459963, 2.69273424634253e-9, 0.927847696699833),
    agn = c(560358.140569773, 1497718.45761725, 174173.819605374,
            15714.4012340279, 7445.49428677381, 2.99792458e-196,
            2.99792458e-198, 2.99792458e-200, 0.92784769400874),
    dustlum = c(3949353348.49075, 15638083468.0756, NA, 19587436816.5664),
    photom = c(7.92153760490619e-7, 0.000199342854726691, 0.00229733486416735,
               0.0203302567233486, 1.18981282095667, 0.138745856682181,
               0.000117685441680569, 3.16100927280049e-7, 545.894815935876)
  ),
  stars = c(114445.181773007, 200401.6216405, 21526.7219048865,
            1604.1112382231, 2690.03349321587, 11.8191908844116,
            0.000100251404459963, 2.69273424634253e-9, 2.69109246507348e-9),
  legacy = c(114445.181773019, 200401.980241726, 21537.8283117996,
             1633.25248356812, 2731.94540274596, 11.8191908844116,
             0.000100251404459963, 2.69273424634253e-9, 1.43292650950679)
)

testthat::test_that('default and explicit viewing-angle output retain the baseline', {
  if (identical(testthat::testing_package(), 'ProSpect')) {
    testthat::expect_identical(.aa_code, asNamespace('ProSpect'))
  } else {
    testthat::expect_identical(isNamespace(.aa_code), FALSE)
  }
  testthat::expect_identical(environment(.aa_code$ProSpectSED), .aa_code)
  for (model in c('Fritz', 'SKIRTOR')) {
    args = .aa_args(model)
    default = .aa_run(args)
    explicit = .aa_run(args, AGN_screen_model = 'viewing_angle')
    testthat::expect_equal(default, explicit, tolerance = 1e-12)
    .aa_finite(default)
    testthat::expect_equal(default$FinalLum$wave, 10^args$waveout)
    testthat::expect_equal(default$FinalLum$lum,
                           .aa_golden[[model]]$final, tolerance = 1e-10)
    testthat::expect_equal(default$AGN$lum,
                           .aa_golden[[model]]$agn, tolerance = 1e-10)
    testthat::expect_equal(unname(default$dustlum),
                           .aa_golden[[model]]$dustlum, tolerance = 1e-10)
    testthat::expect_equal(default$Photom$flux,
                           .aa_golden[[model]]$photom, tolerance = 1e-10)
  }
  args = .aa_args()
  noagn = .aa_run(args, AGN = NULL)
  testthat::expect_equal(noagn$FinalLum$lum, .aa_golden$stars, tolerance = 1e-10)
  legacy = list(Wave = 10^seq(2, 10, length.out = 178),
                Aspec = rep(1e-9, 178))
  old = .aa_run(args, AGN = legacy)
  testthat::expect_equal(old$FinalLum$lum, .aa_golden$legacy, tolerance = 1e-10)
})

testthat::test_that('angle averaging changes reradiation but not the direct spectrum', {
  for (model in c('Fritz', 'SKIRTOR')) {
    args = .aa_args(model)
    mean = .aa_templates(model)$mean
    view = .aa_run(args, AGN_screen_model = 'viewing_angle')
    averaged = .aa_run(args, AGN_screen_model = 'angle_averaged', AGN_anav = mean)
    .aa_finite(averaged)
    testthat::expect_gt(max(abs(
      .aa_on_grid(averaged$AGN$wave, averaged$AGN$lum, view$AGN$wave) -
        view$AGN$lum)), 0)
    testthat::expect_gt(abs(averaged$dustlum['screen'] -
                             view$dustlum['screen']), 0)
    interpolate = .aa_code[[paste0(model, '_anav_interp')]]
    interpolation_args = list(lum = args$AGNlum)
    interpolation_args[paste0(model, '_anav')] = list(mean)
    spectrum = do.call(interpolate, interpolation_args)
    absorbed = sum(spectrum$lum * .aa_code$.erg_to_lsol *
      (1 - exp(-0.3 * (spectrum$wave / 5500)^-0.7)) *
      c(0, diff(spectrum$wave)))
    stars = .aa_stars(tau_birth = 1, tau_screen = 0.3,
                      pow_birth = -0.7, pow_screen = -0.7)
    testthat::expect_equal(unname(averaged$dustlum['screen']),
                           stars$lumtot_screen + absorbed, tolerance = 1e-12)
    testthat::expect_equal(unname(averaged$dustmass['screen']),
                           (stars$lumtot_screen + absorbed) / 2, tolerance = 1e-12)
    interpolate_view = .aa_code[[paste0(model, '_interp')]]
    interpolation_args = list(lum = args$AGNlum)
    interpolation_args[model] = list(args$AGN)
    direct = do.call(interpolate_view, interpolation_args)
    direct_flux = direct$lum * .aa_code$.erg_to_lsol *
      exp(-0.3 * (direct$wave / 5500)^-0.7)
    # Before the radio subtraction (> 1e6 Angstrom), the AGN must be exactly
    # attenuated viewing-angle light plus dust powered by the mean spectrum.
    optical_ir = 10^args$waveout[1:5]
    expected = .aa_on_grid(direct$wave, direct_flux, optical_ir) +
      .aa_on_grid(args$Dale$Wave, args$Dale$Aspec[[1]][1, ] * absorbed, optical_ir)
    testthat::expect_equal(
      .aa_on_grid(averaged$AGN$wave, averaged$AGN$lum, optical_ir),
      expected, tolerance = 1e-12)
    inclined = .aa_run(args, AGN_screen_model = 'angle_averaged',
                       AGN_anav = mean, AGNan = 80)
    testthat::expect_equal(inclined$dustlum, averaged$dustlum, tolerance = 1e-12)
    testthat::expect_gt(abs(inclined$AGN$lum[3] - averaged$AGN$lum[3]), 0)
    zero_view = .aa_run(args, AGN_screen_model = 'viewing_angle', tau_screen = 0)
    zero_mean = .aa_run(args, AGN_screen_model = 'angle_averaged',
                        AGN_anav = mean, tau_screen = 0)
    testthat::expect_equal(zero_mean$FinalLum, zero_view$FinalLum, tolerance = 1e-12)
    testthat::expect_equal(
      .aa_on_grid(zero_mean$AGN$wave, zero_mean$AGN$lum, zero_view$AGN$wave),
      zero_view$AGN$lum, tolerance = 1e-12)
    testthat::expect_equal(unname(zero_mean$dustlum['screen']), 0)
    .aa_finite(zero_mean)
  }
})

testthat::test_that('AGN mode and template classes are validated', {
  args = .aa_args()
  for (mode in list('not-a-mode', NA_character_, character(),
                    c('viewing_angle', 'angle_averaged'), 1)) {
    .aa_error(.aa_run(args, AGN_screen_model = mode), 'arg|screen|mode')
  }
  .aa_error(.aa_run(args, AGN_screen_model = 'angle_averaged',
                    AGN = list(Wave = 1:3, Aspec = 1:3)), 'Fritz|SKIRTOR')
  for (model in c('Fritz', 'SKIRTOR')) {
    args = .aa_args(model)
    other = if (model == 'Fritz') 'SKIRTOR' else 'Fritz'
    .aa_error(.aa_run(args, AGN_screen_model = 'angle_averaged',
                      AGN_anav = .aa_templates(other)$mean), paste0(model, '_anav'))
    .aa_error(.aa_run(args, AGN_screen_model = 'angle_averaged',
                      AGN_anav = unclass(.aa_templates(model)$mean)),
              paste0(model, '_anav'))
    default = .aa_run(args)
    ignored = .aa_run(args, AGN_anav = .aa_templates(other)$mean,
                      AGN_screen_model = 'viewing_angle')
    testthat::expect_equal(default, ignored, tolerance = 1e-12)
  }
  noagn = .aa_run(args, AGNlum = 0, AGN_screen_model = 'angle_averaged')
  testthat::expect_null(noagn$AGN)
  testthat::expect_equal(noagn$FinalLum,
                         .aa_run(args, AGNlum = 0)$FinalLum, tolerance = 1e-12)
})

testthat::test_that('mean interpolators load and validate their own templates', {
  for (model in c('Fritz', 'SKIRTOR')) {
    template_name = paste0(model, '_anav')
    helper_name = paste0(template_name, '_interp')
    mean = .aa_templates(model)$mean
    fixture = .aa_loader_fixture(model)
    helper = fixture$code[[helper_name]]
    supplied = do.call(helper, setNames(list(mean), template_name))
    testthat::expect_equal(length(fixture$calls$names), 0L)
    for (args in list(list(), setNames(list(NULL), template_name))) {
      fixture$calls$names = character()
      fixture$calls$from_helper = logical()
      automatic = do.call(helper, args)
      testthat::expect_identical(fixture$calls$names, template_name, info = model)
      testthat::expect_identical(fixture$calls$from_helper, TRUE, info = model)
      testthat::expect_equal(automatic, supplied, tolerance = 1e-12, info = model)
    }
    other = if (model == 'Fritz') 'SKIRTOR' else 'Fritz'
    for (bad in list(.aa_templates(other)$mean, unclass(mean))) {
      fixture$calls$names = character()
      .aa_error(do.call(helper, setNames(list(bad), template_name)), template_name)
      testthat::expect_equal(length(fixture$calls$names), 0L)
    }
    fixture = .aa_loader_fixture(model, loaded = unclass(mean))
    .aa_error(fixture$code[[helper_name]](), template_name)
    testthat::expect_identical(fixture$calls$names, template_name, info = model)
  }
})

testthat::test_that('ProSpectSED delegates mean loading to the interpolator', {
  for (model in c('Fritz', 'SKIRTOR')) {
    args = .aa_args(model)
    mean = .aa_templates(model)$mean
    fixture = .aa_loader_fixture(model)
    supplied = .aa_run(args, AGN_screen_model = 'angle_averaged',
                       AGN_anav = mean, code = fixture$code)
    testthat::expect_equal(length(fixture$calls$names), 0L)
    automatic = .aa_run(args, AGN_screen_model = 'angle_averaged',
                        code = fixture$code)
    testthat::expect_identical(fixture$calls$names, paste0(model, '_anav'))
    testthat::expect_identical(fixture$calls$from_helper, TRUE)
    testthat::expect_equal(automatic, supplied, tolerance = 1e-12)
    fixture$calls$names = character()
    other = if (model == 'Fritz') 'SKIRTOR' else 'Fritz'
    .aa_error(.aa_run(args, AGN_screen_model = 'angle_averaged',
                      AGN_anav = .aa_templates(other)$mean, code = fixture$code),
              paste0(model, '_anav'))
    for (mode in c('viewing_angle', 'None')) {
      .aa_run(args, AGN_screen_model = mode, code = fixture$code)
    }
    .aa_run(args, AGN_screen_model = 'angle_averaged', AGNlum = 0,
            code = fixture$code)
    testthat::expect_equal(length(fixture$calls$names), 0L)
  }
})

for (.aa_model in c('Fritz', 'SKIRTOR', 'legacy')) {
  for (.aa_grid in c('explicit', 'NULL')) {
    testthat::test_that(paste('None bypasses only the', .aa_model,
                              'ISM screen on the', .aa_grid, 'grid'), {
      model = .aa_model
      args = .aa_none_args(model, .aa_grid)
      for (radio_SF in c(FALSE, TRUE)) {
        for (radio_AGN in c(FALSE, TRUE)) {
          args$addradio_SF = radio_SF
          args$addradio_AGN = radio_AGN
          info = paste('addradio_SF =', radio_SF, 'addradio_AGN =', radio_AGN)
          bypass = .aa_run(args, AGN_screen_model = 'None')
          zero = .aa_run(args, AGN_screen_model = 'viewing_angle', tau_screen = 0)
          screened = .aa_run(args, AGN_screen_model = 'viewing_angle')
          noagn = .aa_run(args, AGNlum = 0)
          .aa_finite(bypass)
          testthat::expect_equal(bypass$AGN, zero$AGN,
                                 tolerance = 1e-12, info = info)
          testthat::expect_gt(max(abs(
            .aa_on_grid(screened$AGN$wave, screened$AGN$lum, bypass$AGN$wave) -
              bypass$AGN$lum)), 0)
          for (name in c('Stars', 'StarsAtten', 'StarsUnAtten', 'DustEmit')) {
            testthat::expect_equal(bypass[[name]], noagn[[name]],
                                   tolerance = 1e-12, info = paste(info, name))
          }
          testthat::expect_gt(max(abs(bypass$StarsAtten$lum -
                                       zero$StarsAtten$lum)), 0)
          testthat::expect_gt(unname(bypass$dustlum['screen']), 0)
          for (name in c('dustlum', 'dustmass')) {
            testthat::expect_equal(bypass[[name]][c('birth', 'screen')],
                                   noagn[[name]][c('birth', 'screen')],
                                   tolerance = 1e-12, info = paste(info, name))
            testthat::expect_gt(unname(screened[[name]]['screen'] -
                                        bypass[[name]]['screen']), 0)
            testthat::expect_equal(bypass[[name]]['AGN'], zero[[name]]['AGN'],
                                   tolerance = 1e-12, info = paste(info, name))
            torus = if (model == 'legacy') unname(zero[[name]]['AGN']) else 0
            if (model == 'legacy') {
              testthat::expect_gt(torus, 0)
            } else {
              testthat::expect_identical(unname(is.na(bypass[[name]]['AGN'])), TRUE)
            }
            testthat::expect_equal(unname(bypass[[name]]['total']),
                                   unname(noagn[[name]]['total']) + torus,
                                   tolerance = 1e-12, info = paste(info, name))
          }
          expected = .aa_on_grid(noagn$FinalLum$wave, noagn$FinalLum$lum,
                                 bypass$FinalLum$wave) +
            .aa_on_grid(zero$AGN$wave, zero$AGN$lum, bypass$FinalLum$wave)
          testthat::expect_equal(bypass$FinalLum$lum, expected,
                                 tolerance = 1e-12, info = info)
          if (model == 'legacy') {
            no_torus = .aa_run(args, AGN_screen_model = 'None', tau_AGN = 0)
            testthat::expect_equal(unname(no_torus$dustlum['AGN']), 0)
            testthat::expect_equal(unname(no_torus$dustmass['AGN']), 0)
            testthat::expect_gt(max(abs(
              .aa_on_grid(no_torus$AGN$wave, no_torus$AGN$lum, bypass$AGN$wave) -
                bypass$AGN$lum)), 0)
          }
        }
      }
    })
  }
}

testthat::test_that('None is exact and remains valid without an active AGN', {
  args = .aa_args()
  for (mode in list('none', 'NONE', NULL, FALSE, c('None', 'viewing_angle'))) {
    .aa_error(.aa_run(args, AGN_screen_model = mode), 'AGN_screen_model.*None')
    .aa_error(.aa_run(args, AGN_screen_model = mode, AGNlum = 0),
              'AGN_screen_model.*None')
  }
  for (inactive in list(list(AGNlum = 0), list(AGN = NULL))) {
    current = args
    current[names(inactive)] = inactive
    testthat::expect_equal(.aa_run(current, AGN_screen_model = 'None'),
                           .aa_run(current), tolerance = 1e-12)
  }
})

testthat::test_that('likelihood forwards top-level and arglist AGN settings', {
  testthat::skip_if_not_installed('ParmOff')
  for (model in c('Fritz', 'SKIRTOR')) {
    args = .aa_args(model)
    mean = .aa_templates(model)$mean
    truth = .aa_run(args, AGN_screen_model = 'angle_averaged', AGN_anav = mean)
    spec = truth$Photom
    spec$fluxerr = pmax(abs(spec$flux) * 0.05, 1e-30)
    fixed = args[setdiff(names(args), c('SFH', 'AGN', 'Dale', 'Dale_M2L_func',
                                      'AGNlum', 'filters', 'returnall'))]
    Data = list(SFH = args$SFH, AGN = args$AGN, Dale = args$Dale,
                Dale_M2L_func = args$Dale_M2L_func, arglist = fixed,
                parm.names = 'AGNlum', logged = FALSE, fit = 'check',
                mode = 'spec', spec = spec, like = 'norm', mon.names = 'LP')
    for (location in c('top-level', 'arglist')) {
      current = Data
      if (location == 'top-level') {
        current$AGN_screen_model = 'angle_averaged'
        current$AGN_anav = mean
        current$arglist$AGN_screen_model = 'viewing_angle'
        other = if (model == 'Fritz') 'SKIRTOR' else 'Fritz'
        current$arglist$AGN_anav = .aa_templates(other)$mean
      } else {
        current$arglist$AGN_screen_model = 'angle_averaged'
        current$arglist$AGN_anav = mean
      }
      checked = .aa_code$ProSpectSEDlike(args$AGNlum, current)
      testthat::expect_equal(checked$SEDout$FinalLum, truth$FinalLum,
                             tolerance = 1e-12, info = location)
      testthat::expect_equal(checked$LP, nrow(spec) * stats::dnorm(0, log = TRUE))
      current$fit = 'optim'
      testthat::expect_equal(.aa_code$ProSpectSEDlike(args$AGNlum, current),
                             -checked$LP, tolerance = 1e-12, info = location)
    }
    invalid = Data
    invalid$AGN_screen_model = 'invalid'
    .aa_error(.aa_code$ProSpectSEDlike(args$AGNlum, invalid), 'screen|arg|mode')
  }
})

testthat::test_that('likelihood forwards None with top-level precedence', {
  testthat::skip_if_not_installed('ParmOff')
  for (model in c('Fritz', 'SKIRTOR', 'legacy')) {
    for (grid in c('explicit', 'NULL')) {
      for (radio in c(FALSE, TRUE)) {
        args = .aa_none_args(model, grid)
        args$addradio_SF = radio
        args$addradio_AGN = radio
        truth = .aa_run(args, AGN_screen_model = 'None')
        spec = truth$Photom
        spec$fluxerr = pmax(abs(spec$flux) * 0.05, 1e-30)
        fixed = args[setdiff(names(args), c('SFH', 'AGN', 'Dale', 'Dale_M2L_func',
                                          'AGNlum', 'filters', 'returnall'))]
        Data = list(SFH = args$SFH, AGN = args$AGN, Dale = args$Dale,
                    Dale_M2L_func = args$Dale_M2L_func, arglist = fixed,
                    parm.names = 'AGNlum', logged = FALSE, fit = 'check',
                    mode = 'spec', spec = spec, like = 'norm', mon.names = 'LP')
        for (location in c('top-level', 'arglist')) {
          current = Data
          if (location == 'top-level') {
            current$AGN_screen_model = 'None'
            current$arglist$AGN_screen_model = 'invalid'
          } else {
            current$arglist$AGN_screen_model = 'None'
          }
          info = paste(model, grid, 'radio =', radio, location)
          checked = .aa_code$ProSpectSEDlike(args$AGNlum, current)
          for (name in c('FinalLum', 'AGN', 'Photom', 'StarsAtten', 'StarsUnAtten',
                         'DustEmit', 'dustlum', 'dustmass')) {
            testthat::expect_equal(checked$SEDout[[name]], truth[[name]],
                                   tolerance = 1e-12, info = paste(info, name))
          }
          testthat::expect_equal(checked$LP,
                                 nrow(spec) * stats::dnorm(0, log = TRUE),
                                 tolerance = 1e-12, info = info)
          current$fit = 'optim'
          testthat::expect_equal(.aa_code$ProSpectSEDlike(args$AGNlum, current),
                                 -checked$LP, tolerance = 1e-12, info = info)
        }
      }
    }
  }
})

testthat::test_that('likelihood validates screen models in both locations', {
  for (location in c('top-level', 'arglist')) {
    for (mode in list('invalid', 'none', NA_character_, character(),
                      c('None', 'viewing_angle'), 1, FALSE)) {
      Data = list(arglist = list())
      if (location == 'top-level') {
        Data$AGN_screen_model = mode
        Data$arglist$AGN_screen_model = 'None'
      } else {
        Data$arglist$AGN_screen_model = mode
      }
      .aa_error(.aa_code$ProSpectSEDlike(1e44, Data), 'AGN_screen_model.*None')
    }
  }
})

for (.aa_model in c('Fritz', 'SKIRTOR')) {
  testthat::test_that(paste('available', .aa_model, 'data exercise real interpolation'), {
    model = .aa_model
    testthat::skip_if_not_installed('ProSpectData')
    available = utils::data(package = 'ProSpectData')$results[, 'Item']
    datasets = c(model, paste0(model, '_anav'))
    if (!all(datasets %in% available)) {
      testthat::skip(paste('ProSpectData is missing', paste(datasets, collapse = '/')))
    }
    env = new.env(parent = emptyenv())
    for (name in datasets) {
      utils::data(list = name, package = 'ProSpectData', envir = env)
    }
    args = .aa_args(model)
    args$AGN = env[[model]]
    view = .aa_run(args)
    averaged = .aa_run(args, AGN_screen_model = 'angle_averaged',
                       AGN_anav = env[[paste0(model, '_anav')]])
    automatic = .aa_run(args, AGN_screen_model = 'angle_averaged')
    testthat::expect_equal(automatic$FinalLum, averaged$FinalLum, tolerance = 1e-12)
    .aa_finite(view)
    .aa_finite(averaged)
    testthat::expect_gt(max(abs(averaged$FinalLum$lum - view$FinalLum$lum)), 0)
    zero_view = .aa_run(args, tau_screen = 0)
    zero_mean = .aa_run(args, tau_screen = 0, AGN_screen_model = 'angle_averaged',
                        AGN_anav = env[[paste0(model, '_anav')]])
    testthat::expect_equal(zero_mean$FinalLum, zero_view$FinalLum, tolerance = 1e-10)
  })
}
