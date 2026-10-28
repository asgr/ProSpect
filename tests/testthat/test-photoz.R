# Photo-z fitting mode of ProSpectSEDlike().
#
# The branch that builds the cosmology used to read Data$arglist$HO rather than
# $H0, so a photo-z fit that supplied H0/OmegaM/OmegaL (rather than a named
# 'ref') silently passed H0 = NULL to celestial and produced an unusable model.
# These tests drive the branch directly.

.ps_photoz_data = function(truth_z = 0.3, cosmology = NULL, ref = NULL,
                           z_genSF = NULL) {
  bands = c('g_SDSS', 'r_SDSS', 'K_VISTA')
  filtout = .ps_filtout(bands)

  # Generate the mock at the truth redshift, using the same cosmology.
  args = .ps_args(filtout = filtout, z = truth_z, sparse = 10)
  if (!is.null(ref)) args$ref = ref
  if (!is.null(cosmology)) args = c(args, cosmology)
  flux = do.call(ProSpectSED, args)$Photom

  arglist = list(photoz = TRUE, massfunc = massfunc_b5, IGMfunc = NULL)
  if (!is.null(cosmology)) arglist = c(arglist, cosmology)
  if (!is.null(ref)) arglist$ref = ref
  if (!is.null(z_genSF)) arglist$z_genSF = z_genSF

  list(
    truth_z = truth_z,
    Data = list(
      flux = data.frame(filter = bands, pivwave = 1, flux = flux,
                        fluxerr = flux * 0.05),
      filtout = filtout,
      parm.names = 'z',
      speclib = .ps_speclib(),
      Dale = .ps_dale(),
      SFH = SFHfunc,
      fit = 'optim',
      arglist = arglist
    )
  )
}

test_that('photo-z fitting works when H0 is supplied without a ref', {
  skip_if_not_installed('celestial')

  mock = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692))

  # This used to error because H0 was read from the non-existent arglist$HO.
  value = ProSpectSEDlike(mock$truth_z, mock$Data)
  expect_type(value, 'double')
  expect_true(is.finite(value))
})

test_that('photo-z likelihood is minimised at the true redshift', {
  skip_if_not_installed('celestial')

  mock = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692))
  truth = mock$truth_z

  at_truth = ProSpectSEDlike(truth, mock$Data)
  perturbed = vapply(c(truth - 0.1, truth + 0.1, truth + 0.3),
                     function(z) ProSpectSEDlike(z, mock$Data), numeric(1))

  expect_true(all(perturbed > at_truth))
})

test_that('photo-z fitting gives the same answer via ref or explicit cosmology', {
  skip_if_not_installed('celestial')

  ref_mock = .ps_photoz_data(ref = 'Planck18')
  cosmo_mock = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692))

  # Planck18 and these parameters are nearly the same cosmology, so the
  # likelihoods should agree closely (not necessarily exactly).
  expect_equal(
    ProSpectSEDlike(ref_mock$truth_z, ref_mock$Data),
    ProSpectSEDlike(cosmo_mock$truth_z, cosmo_mock$Data),
    tolerance = 0.05
  )
})

test_that('photo-z replaces agemax and LumDist_Mpc from the fitted redshift', {
  skip_if_not_installed('celestial')

  mock = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692))
  checked = ProSpectSEDlike(mock$truth_z, utils::modifyList(mock$Data, list(fit = 'check')))

  # The branch must have written a finite, positive distance into the arglist.
  expect_true(is.finite(checked$Data$arglist$LumDist_Mpc))
  expect_gt(checked$Data$arglist$LumDist_Mpc, 0)

  # And it must match the celestial value for the fitted redshift.
  expect_equal(
    checked$Data$arglist$LumDist_Mpc,
    celestial::cosdistLumDist(mock$truth_z, H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692),
    tolerance = 1e-8
  )
})

test_that('photo-z with z_genSF shortens agemax by the age at formation', {
  skip_if_not_installed('celestial')

  plain = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692))
  # z_genSF must exceed the fitted redshift: stars cannot start forming after
  # the epoch being observed.
  delayed = .ps_photoz_data(cosmology = list(H0 = 67.8, OmegaM = 0.308, OmegaL = 0.692),
                            z_genSF = 2)

  plain_check = ProSpectSEDlike(plain$truth_z,
                                utils::modifyList(plain$Data, list(fit = 'check')))
  delayed_check = ProSpectSEDlike(delayed$truth_z,
                                  utils::modifyList(delayed$Data, list(fit = 'check')))

  # Star formation starting at z_genSF = 0.1 leaves less time before z = 0.3
  # than starting at the Big Bang, so agemax must be smaller.
  expect_lt(delayed_check$Data$arglist$agemax, plain_check$Data$arglist$agemax)
  expect_gt(delayed_check$Data$arglist$agemax, 0)
})

test_that('photo-z demands a cosmology when no ref is given', {
  skip_if_not_installed('celestial')

  mock = .ps_photoz_data()
  mock$Data$arglist$OmegaM = NULL

  expect_error(ProSpectSEDlike(mock$truth_z, mock$Data),
               'Cosmology must be specified for photo-z')
})