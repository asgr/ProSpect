# Filter-registry integrity. getfilt() is the entry point for every photometric
# calculation, so a missing or malformed filter silently poisons photometry.
# Filter curves are two-column tables in data/: wave in Angstrom, response.

test_that('getfilt returns a valid wave/response table for common bands', {
  bands = c('FUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'K_VISTA', 'W1_WISE',
            'P100_Herschel', 'S500_Herschel')

  for (band in bands) {
    curve = getfilt(band)

    expect_false(is.null(curve), label = paste(band, 'resolves'))
    expect_equal(ncol(curve), 2)
    expect_equal(nrow(curve) >= 2, TRUE)
    # Wavelengths ascend, and responses are non-negative.
    expect_true(all(diff(curve[, 1]) > 0), label = paste(band, 'ascending wave'))
    expect_true(all(curve[, 2] >= 0), label = paste(band, 'non-negative response'))
    # The curve must actually respond to something.
    expect_gt(max(curve[, 2]), 0, label = paste(band, 'has a response'))
  }
})

test_that('getfilt accepts the documented shorthands', {
  # Each shorthand must resolve to the same curve as the full name.
  pairs = list(
    c('W1', 'W1_WISE'),
    c('I1', 'I1_Spitzer'),
    c('P100', 'P100_Herschel'),
    c('S850', 'S850_JCMT'),
    c('Band3', 'Band3_ALMA'),
    c('BandC', 'BandC_VLA')
  )

  for (pair in pairs) {
    expect_equal(getfilt(pair[1]), getfilt(pair[2]),
                 label = paste(pair[1], 'resolves to', pair[2]))
  }
})

test_that('getfilt builds tophat filters on demand', {
  top = getfilt('tophat_4000_6000')

  expect_equal(top[, 1], c(4000, 6000))
  expect_equal(top[, 2], c(1, 1))
})

test_that('getfilt messages rather than errors on an unknown filter', {
  # A typo returns NULL with a message. Callers must therefore check for NULL;
  # this test documents that contract rather than endorsing it.
  expect_message(out <- getfilt('not_a_real_filter'), 'not recognised')
  expect_null(out)
})

# Pull the valid_filters vector out of getfilt()'s body. The registry is a
# local variable rather than an exported object, and walking the AST is far
# more robust than deparsing and regexing the source.
.valid_filters = function() {
  found = list()
  walk = function(node) {
    if (is.call(node)) {
      if (as.character(node[[1]]) %in% c('=', '<-') &&
          identical(as.character(node[[2]]), 'valid_filters')) {
        found[[1]] <<- node[[3]]
        return(invisible(NULL))
      }
      for (i in seq_along(node)[-1]) walk(node[[i]])
    }
    invisible(NULL)
  }
  walk(body(getfilt))
  if (length(found) == 0) return(character(0))
  eval(found[[1]])
}

test_that('every filter named in valid_filters resolves', {
  registry = .valid_filters()
  expect_gt(length(registry), 100)

  # `missing` means "cannot be turned into a filter curve", whether that shows
  # up as NULL (installed package) or an error from get() after data() fails
  # (devtools::load_all).
  missing = vapply(registry, function(band) {
    out = tryCatch(
      suppressWarnings(suppressMessages(getfilt(band))),
      error = function(e) NULL
    )
    is.null(out)
  }, logical(1))

  expect_equal(registry[missing], character(0))
})

# All \alias{} entries from the package documentation, as a named character
# vector mapping alias -> Rd file. Read from the built help index: under
# R CMD check the package is installed and man/*.Rd is not on disk, so parsing
# the source tree would fail there. tools::Rd_db() needs the Rd sources too,
# whereas help/aliases.rds is present in every installed package.
.doc_aliases = function() {
  path = system.file('help', 'aliases.rds', package = 'ProSpect')
  if (!nzchar(path)) {
    # devtools::load_all() on a source tree that has never been installed
    rd = list.files(testthat::test_path('..', '..', 'man'), pattern = '\\.Rd$',
                    full.names = TRUE)
    return(unlist(lapply(rd, function(file) {
      sub('^\\\\alias\\{(.*)\\}$', '\\1',
          grep('^\\\\alias\\{', readLines(file, warn = FALSE), value = TRUE))
    })))
  }
  names(readRDS(path))
}

# Names of the filter curves shipped with the package. Query the data index
# rather than the data/ directory: with LazyData the files are compiled into
# Rdata.rdb and there is no filt_*.tab on disk in the installed package.
.filt_data_items = function() {
  items = utils::data(package = 'ProSpect')$results[, 'Item']
  sub('^filt_', '', grep('^filt_', items, value = TRUE))
}

test_that('every ALMA band resolves, including Band 1', {
  # Band1_ALMA used to be registered as " Band1_ALMA" with a leading space, so
  # both the shorthand and the full name failed even though the data object
  # existed. Check the whole family now works, by shorthand and full name.
  bands = paste0('Band', 1:9)

  for (band in bands) {
    curve = suppressWarnings(suppressMessages(getfilt(band)))
    expect_false(is.null(curve), label = paste(band, 'resolves'))
    expect_equal(curve, suppressWarnings(suppressMessages(getfilt(paste0(band, '_ALMA')))),
                 label = paste(band, 'shorthand matches full name'))
  }

  # Band 1 is the longest wavelength of the ALMA set (3 mm).
  peaks = vapply(bands, function(band) {
    curve = suppressWarnings(suppressMessages(getfilt(band)))
    curve[which.max(curve[, 2]), 1]
  }, numeric(1))
  expect_equal(names(peaks)[which.max(peaks)], 'Band1')
})

test_that('every filt_* data object is registered in getfilt', {
  # Catches the reverse mistake: adding a filter curve but forgetting to list it
  # in valid_filters.
  bands = .filt_data_items()
  expect_gt(length(bands), 50)

  unresolved = vapply(bands, function(band) {
    is.null(tryCatch(suppressWarnings(suppressMessages(getfilt(band))),
                     error = function(e) NULL))
  }, logical(1))

  expect_equal(bands[unresolved], character(0))
})

test_that('Rubin filters alias the LSST curves', {
  # Rubin Observatory is the renamed LSST, so there is one curve per band and
  # the Rubin name is an alias rather than a separate data object.
  for (band in c('u', 'g', 'r', 'i', 'z', 'y')) {
    rubin = suppressWarnings(suppressMessages(getfilt(paste0(band, '_Rubin'))))
    lsst = suppressWarnings(suppressMessages(getfilt(paste0(band, '_LSST'))))

    expect_false(is.null(rubin), label = paste(band, 'Rubin resolves'))
    expect_equal(rubin, lsst, label = paste(band, 'Rubin matches LSST'))
  }

  # No separate Rubin data objects are shipped, so this guards against someone
  # adding duplicate curves that would then drift apart.
  expect_equal(grep('_Rubin$', .filt_data_items(), value = TRUE), character(0))
})

test_that('the y band is resolvable under either capitalisation per survey', {
  # Each survey has a single y curve, so y_/Y_ must not be distinct filters.
  expect_equal(suppressMessages(getfilt('y_LSST')), suppressMessages(getfilt('Y_LSST')))
  expect_equal(suppressMessages(getfilt('y_HSC')), suppressMessages(getfilt('Y_HSC')))
  expect_equal(suppressMessages(getfilt('y_Rubin')), suppressMessages(getfilt('y_LSST')))
  expect_equal(suppressMessages(getfilt('Y_Rubin')), suppressMessages(getfilt('y_LSST')))

  # LSST y and HSC Y are genuinely different curves, so they must not alias.
  expect_false(identical(suppressMessages(getfilt('y_LSST')),
                         suppressMessages(getfilt('Y_HSC'))))
})

test_that('pivwave and cenwave cover the broadband filters', {
  e = new.env()
  data('pivwave', envir = e)
  data('cenwave', envir = e)

  expect_named(e$pivwave, c('filter', 'pivwave'))
  expect_named(e$cenwave, c('filter', 'cenwave'))
  expect_equal(nrow(e$pivwave), nrow(e$cenwave))

  # The optical/NIR/FIR broadbands used throughout the vignettes must be there,
  # since filters = 'all' is built from cenwave.
  for (band in c('FUV_GALEX', 'u_SDSS', 'K_VISTA', 'W1_WISE', 'P100_Herschel',
                 'S500_Herschel')) {
    expect_true(band %in% e$pivwave$filter, label = paste(band, 'in pivwave'))
    expect_true(band %in% e$cenwave$filter, label = paste(band, 'in cenwave'))
  }

  # Pivot and central wavelengths must be physical and sensibly ordered.
  expect_true(all(e$pivwave$pivwave > 0))
  expect_true(all(e$cenwave$cenwave > 0))
})

test_that('filters = "all" resolves to a usable filter list', {
  # ProSpectSED builds filtout from cenwave when filters = 'all'.
  sed = ProSpectSED(massfunc = massfunc_b5, speclib = .ps_speclib(),
                    Dale = .ps_dale(), filters = 'all', z = 0.1, sparse = 10)

  expect_true(length(sed$filtout) > 50)
  expect_false(any(vapply(sed$filtout, is.null, logical(1))))
  expect_equal(length(sed$Photom), length(sed$filtout))
  expect_true(all(is.finite(sed$Photom)))
})
test_that('every filter data object is documented with an alias', {
  # R CMD check warns twice when a data object has no \alias{}: once for missing
  # documentation and once for an undocumented data set. The nine Roman WFI
  # filters used to fail both because man/filters.Rd said filt_*_ROMAN while the
  # objects are filt_*_WFI_ROMAN.
  aliases = .doc_aliases()
  expect_gt(length(aliases), 100)

  undocumented = setdiff(paste0('filt_', .filt_data_items()), aliases)
  expect_equal(undocumented, character(0))
})

test_that('no documentation aliases a filter object that does not exist', {
  # The mirror of the test above: a stale \alias{} makes R CMD check complain
  # that data are used in an Rd file but are not in the code.
  aliases = grep('^filt_', .doc_aliases(), value = TRUE)
  missing = setdiff(aliases, paste0('filt_', .filt_data_items()))

  expect_equal(missing, character(0))
})

test_that('the Roman WFI filters resolve under their data object names', {
  # Nine current-generation Roman filters, distinct from the older WFIRST set.
  wfi = c('F062', 'F087', 'F106', 'F129', 'F146', 'F158', 'F184', 'F213', 'W146')

  for (band in wfi) {
    name = paste0(band, '_WFI_ROMAN')
    curve = suppressWarnings(suppressMessages(getfilt(name)))

    expect_false(is.null(curve), label = paste(name, 'resolves'))
    expect_equal(ncol(curve), 2)
    expect_true(all(diff(curve[, 1]) > 0), label = paste(name, 'ascending wave'))
  }

  # The two sets use different naming conventions, which is why they look
  # unrelated: WFIRST borrowed the HST-style R062/Z087/Y106/J129/H158 names,
  # while the WFI set uses the Roman filter names F062/F087/F106/F129/F158.
  # Only F184 and W146 are named identically in both.
  wfirst = grep('_WFIRST$', .filt_data_items(), value = TRUE)
  expect_length(wfirst, 7)
  wfirst_bands = sub('_WFIRST$', '', wfirst)
  expect_setequal(intersect(wfi, wfirst_bands), c('F184', 'W146'))
})

test_that('the Roman WFI and WFIRST filter sets are genuinely different curves', {
  # They share several band names but are NOT aliases: different response
  # curves over different wavelength ranges. Documentation says so, and this
  # stops someone "helpfully" aliasing them together later.
  shared = intersect(sub('_WFIRST$', '', grep('_WFIRST$', .filt_data_items(), value = TRUE)),
                     sub('_WFI_ROMAN$', '', grep('_WFI_ROMAN$', .filt_data_items(), value = TRUE)))
  expect_setequal(shared, c('F184', 'W146'))

  for (band in shared) {
    wfirst = suppressWarnings(suppressMessages(getfilt(paste0(band, '_WFIRST'))))
    wfi = suppressWarnings(suppressMessages(getfilt(paste0(band, '_WFI_ROMAN'))))

    expect_false(identical(wfirst, wfi), label = paste(band, 'sets differ'))
  }
})
