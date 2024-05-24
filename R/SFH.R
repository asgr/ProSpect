SFHfunc = function(massfunc = massfunc_b5,
                   forcemass = FALSE,
                   agescale = 1,
                   stellpop = 'BC03lr',
                   speclib = NULL,
                   tau_birth = 1.0,
                   tau_screen = 0.3,
                   pow_birth = -0.7,
                   pow_screen = -0.7,
                   filters = 'all',
                   Z = 5,
                   emission = FALSE,
                   veldisp = 50,
                   range = 5,
                   res = 0.5,
                   emission_scale = 'FUV',
                   escape_frac = 1 - emission,
                   Ly_limit = 911.75,
                   LKL10 = NULL,
                   disp_stars = FALSE,
                   LSF = NULL,
                   z = 0.1,
                   H0 = 67.8,
                   OmegaM = 0.308,
                   OmegaL = 1 - OmegaM,
                   ref,
                   outtype = 'mag',
                   sparse = 5,
                   intSFR = FALSE,
                   unimax = 13.8e9,
                   agemax = NULL,
                   LumDist_Mpc = NULL,
                   Eb = 0,
                   L0 = 2175.8,
                   LFWHM = 470,
                   ...) {
  #Ly_limit should be 911.75 Ang (the actual ionisation limit) or sometimes 1215.67 Ang (Lyman alpha)

  dots = list(...)
  massfunc_args = dots[names(dots) %in% names(formals(massfunc))]

  if (is.null(speclib)) {
    if (stellpop == 'BC03lr') {
      BC03lr = NULL
      data('BC03lr', envir = environment())
      speclib = BC03lr
    }else if (stellpop == 'BC03hr') {
      BC03hr = NULL
      data('BC03hr', envir = environment())
      speclib = BC03hr
    }else if (stellpop == 'EMILES') {
      EMILES = NULL
      data('EMILES', envir = environment())
      speclib = EMILES
    }else if (stellpop == 'BPASS') {
      BPASS = NULL
      data('BPASS', envir = environment())
      speclib = BPASS
    }else{
      stop('Need speclib or stellpop!')
    }
  }

  if (any(speclib$Age <= 1e7)) {
    birthcloud = max(which(speclib$Age <= 1e7))
  } else{
    birthcloud = 1
  }

  if (!is.null(filters)) {
    if (filters[1] == 'all') {
      cenwave = NULL
      data('cenwave', envir = environment())
      filters = cenwave$filter
    }else if(filters[1] == 'GAMA'){
      filters = c(
        'FUV_GALEX',
        'NUV_GALEX',
        'u_VST',
        'g_VST',
        'r_VST',
        'i_VST',
        'Z_VISTA',
        'Y_VISTA',
        'J_VISTA',
        'H_VISTA',
        'K_VISTA',
        'W1_WISE' ,
        'W2_WISE',
        'W3_WISE',
        'W4_WISE',
        'P100_Herschel',
        'P160_Herschel',
        'S250_Herschel' ,
        'S350_Herschel',
        'S500_Herschel'
      )
    }else if(filters[1] == 'WAVES'){
      filters = c(
        'u_VST',
        'g_VST',
        'r_VST',
        'i_VST',
        'Z_VISTA',
        'Y_VISTA',
        'J_VISTA',
        'H_VISTA',
        'K_VISTA',
        'W1_WISE' ,
        'W2_WISE'
      )
    }
  }

  if (!is.function(Z)) {
    if (Z %% 1 != 0) {
      #Check if the Z is non integer, if so then convert to a function
      tempZfunc = function(age, Z, ...) {
        rep(Z, length(age))
      }
      formals(tempZfunc)$Z = Z
      Z = tempZfunc
    }
  }

  # if(requireNamespace('Rfast')){
  #   #this proved to be slower in fact, but leave commented for future tests
  #   colsum_loc = Rfast::colsums
  #   formals(colsum_loc)$parallel=TRUE
  # }else{
  #   colsum_loc = colSums
  # }

  if(agescale != 1){
    agevec = speclib$Age * agescale
  }else{
    agevec = speclib$Age
  }

  if (is.function(Z)) {
    dots = list(...)
    Z_args = dots[names(dots) %in% names(formals(Z))]
    Zvec = do.call('Z', c(
      list(agevec),
      list(massfunc = massfunc),
      Z_args,
      massfunc_args
    ))
    Zlist = interp_param(Zvec, speclib$Z, log = TRUE)
    Zwmat = matrix(0, length(speclib$Age), length(speclib$Z))
    Zwmat[cbind(1:length(speclib$Age), Zlist[, 'ID_hi'])] = Zlist[, 'wt_hi']
    Zwmat[cbind(1:length(speclib$Age), Zlist[, 'ID_lo'])] = Zlist[, 'wt_lo']
    Zuse = which(colSums(Zwmat) > 0)
    Zdoweight = TRUE
  } else{
    Zuse = Z
    Zvec = rep(speclib$Z[Zuse], length(speclib$Age))
    Zdoweight = FALSE
  }

  wave_lum = speclib$Wave
  if (sparse > 1) {
    sparse = seq(1, dim(speclib$Zspec[[1]])[2], by = sparse)
    for (i in Zuse) {
      speclib$Zspec[[i]] = speclib$Zspec[[i]][, sparse]
    }
    wave_lum = wave_lum[sparse]
  }

  if (intSFR) {
    massvec = rep(0, length(agevec))
    for (i in 1:length(agevec)) {
      tempint = try(do.call('integrate', c(
        list(
          f = massfunc,
          lower = speclib$AgeBins[i] * agescale,
          upper = speclib$AgeBins[i + 1] * agescale
        ),
        massfunc_args
      ))$value,
      silent = TRUE)
      if (inherits(tempint, "try-error")) {
        massvec[i] = 0
      } else{
        massvec[i] = tempint
      }
    }
    if (sum(massvec, na.rm = TRUE) == 0) {
      #pracma integral seems to work for very bursty star formation where integrate fails
      massvec = rep(0, length(agevec))
      for (i in 1:length(agevec)) {
        tempint = try(do.call('integral', c(
          list(
            f = massfunc,
            xmin = speclib$AgeBins[i] * agescale,
            xmax = speclib$AgeBins[i + 1] * agescale
          ),
          massfunc_args
        )),
        silent = TRUE)
        if (inherits(tempint, "try-error")) {
          massvec[i] = 0
        } else{
          massvec[i] = tempint
        }
      }
    }
  } else{
    massvec = do.call('massfunc', c(list(agevec), massfunc_args)) * speclib$AgeWeights
  }

  if (unimax != FALSE) {
    if (is.null(agemax)) {
      agemax = unimax - cosdistTravelTime(
        z = z,
        H0 = H0,
        OmegaM = OmegaM,
        OmegaL = OmegaL,
        ref = ref
      ) * 1e9
    }
  }

  if (!is.null(agemax)) {
    massvec[speclib$Age > agemax] = 0
  }

  if (any(massvec < 0)) {
    stop('Supplied massfunc cannot create negative SFR!')
  }

  if (forcemass == FALSE) {
    forcescale = 1
    masstot = sum(massvec)
  } else{
    masstot = sum(massvec)
    forcescale = forcemass / masstot
    massvec = massvec * forcescale
    masstot = forcemass
  }

  # Here we modify the speclib if we have an escape fraction less than 1.
  # This is because this will be the first process that occurs, before birth cloud dust or ISM screen dust
  # This means we also need to track the luminosity of the stars before any UV absorption (lum)

  if (length(Zuse) > 1) {
    lum = rep(0, length(wave_lum))
    for (Zid in Zuse) {
      toadd = colSums(speclib$Zspec[[Zid]] * massvec * Zwmat[, Zid])
      #toadd[which(is.na(toadd) | (toadd < 0))] = 0 #since some spectral libraries have negative flux (noisy empirical data)
      lum = lum + toadd
      if (any(escape_frac < 1)) {
        if (length(Ly_limit) == 1) {
          speclib$Zspec[[Zid]][, wave_lum < Ly_limit] = speclib$Zspec[[Zid]][, wave_lum < Ly_limit] * escape_frac
        } else{
          for (i in 1:(length(Ly_limit) - 1)) {
            sel = which(wave_lum < Ly_limit[i] & wave_lum > Ly_limit[i + 1])
            speclib$Zspec[[Zid]][, sel] = speclib$Zspec[[Zid]][, sel] *
              escape_frac[i]
          }
          sel = which(wave_lum < Ly_limit[length(Ly_limit)])
          speclib$Zspec[[Zid]][, sel] = speclib$Zspec[[Zid]][, sel] * escape_frac[length(Ly_limit)]
        }
      }
    }
  } else{
    lum = colSums(speclib$Zspec[[Zuse]] * massvec)
    # if(any(escape_frac<1)){
    #   speclib$Zspec[[Zuse]][,wave_lum<Ly_limit]=speclib$Zspec[[Zuse]][,wave_lum<Ly_limit]*escape_frac
    # }
    if (any(escape_frac < 1)) {
      if (length(Ly_limit) == 1) {
        speclib$Zspec[[Zuse]][, wave_lum < Ly_limit] = speclib$Zspec[[Zuse]][, wave_lum < Ly_limit] * escape_frac
      } else{
        for (i in 1:(length(Ly_limit) - 1)) {
          sel = which(wave_lum < Ly_limit[i] & wave_lum > Ly_limit[i + 1])
          speclib$Zspec[[Zuse]][, sel] = speclib$Zspec[[Zuse]][, sel] *
            escape_frac[i]
        }
        sel = which(wave_lum < Ly_limit[length(Ly_limit)])
        speclib$Zspec[[Zuse]][, sel] = speclib$Zspec[[Zuse]][, sel] * escape_frac[length(Ly_limit)]
      }
    }
  }

  # This should be pre dispersion being added, and birthcloud or ISM attenuation.
  lumtot_unatten = sum(.qdiff(wave_lum) * lum)
  lum_unatten = lum

  if (tau_birth != 0) {
    lum = rep(0, length(wave_lum))
    for (Zid in Zuse) {
      if (tau_birth != 0) {
        #speclib$Zspec[[Zid]][1:birthcloud,]=t(t(speclib$Zspec[[Zid]][1:birthcloud,])*CF_birth(wave_lum, tau=tau_birth, pow=pow_birth))
        speclib$Zspec[[Zid]][1:birthcloud, ] = speclib$Zspec[[Zid]][1:birthcloud, ] *
          rep(CF_birth(wave_lum, tau = tau_birth, pow = pow_birth),
              each = birthcloud)
      }
      if (Zdoweight) {
        lum = lum + colSums(speclib$Zspec[[Zid]] * massvec * Zwmat[, Zid])
      } else{
        lum = colSums(speclib$Zspec[[Zid]] * massvec)
      }
    }
    lumtot_birth = lumtot_unatten - sum(.qdiff(wave_lum) * lum)
  }else{
    lumtot_birth = 0
  }

  if(disp_stars){
    # Here we optionally disperse the spectrum along our LoS.
    # Remember here we are still in instrinsic luminosity space

    grid = seq(-range, range, by=res)
    weights = dnorm(grid)
    wave_lum_log = log10(wave_lum)

    if(!is.null(LSF)){
      if(is.function(LSF)){
        vel_LSF = LSF(wave_lum*(1 + z)) #to get LSF dispersion in km/s into z in the oberved frame
      }else if(is.matrix(LSF) | is.data.frame(LSF)){
        vel_LSF = approx(x=log10(LSF[,1]), y=LSF[,2], xout=log10(wave_lum*(1 + z)), rule=2)$y
      }else if(length(LSF == 1)){
        vel_LSF = rep(LSF, length(wave_lum))
      }else{
        stop('LSF is in the wrong format!')
      }
    }else{
      vel_LSF = rep(0, length(wave_lum))
    }

    z_disp = sqrt(veldisp^2 + vel_LSF^2)/(.c_to_mps/1000) #this will be a vector of length wave_lum

    lum_conv = numeric(length(lum))

    for(i in seq_along(grid)){
      z_seq = grid[i]*z_disp
      new_wave = log10(wave_lum*(1 + z_seq))
      new_lum = log10(lum/(1 + z_seq))
      new_lum = 10^approx(x=new_wave, y=new_lum, xout=wave_lum_log, rule=2, yleft=new_lum[1], yright=new_lum[length(new_lum)])$y
      new_lum = new_lum*weights[i]
      lum_conv = lum_conv + new_lum
    }

    lum = lum_conv*res
    rm(lum_conv)
  }

  if (emission) {
    if (emission_scale == 'FUV') {
      if (length(Ly_limit) == 1 | all(escape_frac == escape_frac[1])) {
        sel = which(wave_lum < Ly_limit[1])
        All_lum = (1 - escape_frac) * sum(.qdiff(wave_lum[sel]) * lum_unatten[sel])
      } else{
        All_lum = 0
        for (i in 1:(length(Ly_limit) - 1)) {
          sel = which(wave_lum < Ly_limit[i] & wave_lum > Ly_limit[i + 1])
          All_lum = All_lum + (1 - escape_frac[i]) * (Ly_limit[i] - Ly_limit[i + 1]) * mean(lum_unatten[sel])
        }
        sel = which(wave_lum < Ly_limit[length(Ly_limit)])
        All_lum = All_lum + (1 - escape_frac[length(Ly_limit)]) * sum(.qdiff(wave_lum[sel]) * lum_unatten[sel])
      }

      emission_input = list(All_lum = All_lum,
                            veldisp = veldisp,
                            Z = Zvec[1])
      emissionadd_unatten = emissionLines(All_lum = All_lum,
                                          veldisp = veldisp,
                                          Z = Zvec[1],
                                          LSF = LSF,
                                          z_LSF = z,
                                          LKL10 = LKL10,
                                          res = res,
                                          range = range)
    } else if (emission_scale == 'SFR') {
      SFRburst_emission = (1 - escape_frac) * do.call('integrate', c(list(
        f = massfunc,
        lower = 0,
        upper = 1e7
      ), massfunc_args))$value * forcescale / 1e7
      emission_input = list(SFR = SFRburst_emission,
                            veldisp = veldisp,
                            Z = Zvec[1])
      emissionadd_unatten = emissionLines(SFR = SFRburst_emission,
                                          veldisp = veldisp,
                                          Z = Zvec[1],
                                          LSF = LSF,
                                          z_LSF = z,
                                          LKL10 = LKL10,
                                          res = res,
                                          range = range)
    } else{
      stop('emission_scale must be one of SFR or FUV!')
    }

    emissionadd_atten = emissionadd_unatten
    if(tau_birth != 0){
      emissionadd_atten$lum = emissionadd_atten$lum * CF_birth(emissionadd_atten$wave, tau = tau_birth, pow = pow_birth)
    }

    lumtot_emission_unatten = sum(.qdiff(emissionadd_unatten$wave) * emissionadd_unatten$lum)
    lumtot_emission_atten = sum(.qdiff(emissionadd_atten$wave) * emissionadd_atten$lum)

    lumtot_birth = lumtot_birth - lumtot_emission_unatten + lumtot_emission_atten

    if(lumtot_birth < 0){
     lumtot_birth = 0
    }

    lum = addspec(wave_lum,
                  lum,
                  emissionadd_atten$wave,
                  emissionadd_atten$lum)
    lum_unatten = approxfun(log10(wave_lum), lum_unatten)(log10(lum$wave))

    wave_lum = lum$wave
    lum = lum$flux
    #wave_lum=lum_unatten$wave
    #lum_unatten=lum_unatten$flux
  } else{
    emission_input = NULL
  }

  if (tau_screen != 0) {
    lum = lum * CF_screen(
      wave_lum,
      tau = tau_screen,
      pow = pow_screen,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM
    )
    lumtot_screen = (lumtot_unatten - lumtot_birth) - sum(.qdiff(wave_lum) * lum)

    if(lumtot_screen < 0){
      lumtot_screen = 0
    }

  } else{
    lumtot_screen = 0
  }

  SFRburst = do.call('integrate', c(list(
    f = massfunc, lower = 0, upper = 1e8
  ), massfunc_args))$value * forcescale / 1e8

  lumtot_atten = lumtot_unatten - sum(.qdiff(wave_lum) * lum)

  if (z < 0 | is.null(filters)) {
    out = NULL
    flux = NULL
  } else{
    if (z > 0) {
      flux = Lum2Flux(
        wave = wave_lum,
        lum = lum,
        z = z,
        H0 = H0,
        OmegaM = OmegaM,
        OmegaL = OmegaL,
        ref = ref,
        LumDist_Mpc = LumDist_Mpc
      )
      if (!is.null(outtype)) {
        out = photom_flux(flux, outtype = outtype, filters = filters)
        if (is.list(filters)) {
          cenout = {
          }
          for (i in filters) {
            cenout = c(cenout, cenwavefunc(i))
          }
          if (all(!is.null(names(filters)))) {
            out = data.frame(
              filter = names(filters),
              cenwave = cenout,
              out = out
            )
          } else{
            out = data.frame(filter = NA,
                             cenwave = cenout,
                             out = out)
          }
        } else{
          cenwave = NULL
          data('cenwave', envir = environment())
          out = data.frame(cenwave[match(filters, cenwave$filter), ], out = out)
        }
      } else{
        out = NULL
      }
    } else{
      flux = cbind(wave = wave_lum, flux = lum * .lsol_to_absolute)
      if (!is.null(outtype)) {
        out = photom_flux(flux, outtype = outtype, filters = filters)
        if (is.list(filters)) {
          cenout = {
          }
          for (i in filters) {
            cenout = c(cenout, cenwavefunc(i))
          }
          if (all(!is.null(names(filters)))) {
            out = data.frame(
              filter = names(filters),
              cenwave = cenout,
              out = out
            )
          } else{
            out = data.frame(filter = NA,
                             cenwave = cenout,
                             out = out)
          }
        } else{
          cenwave = NULL
          data('cenwave', envir = environment())
          out = data.frame(cenwave[match(filters, cenwave$filter), ], out = out)
        }
      } else{
        out = NULL
      }
    }
  }

  SFR = massvec / speclib$AgeWeights

  if (!is.null(agemax)) {
    sel = which(speclib$Age <= agemax)
    agevec = agevec[sel]
    massvec = massvec[sel]
    SFR = SFR[sel]
    Zvec = Zvec[sel]
  }

  return(
    list(
      flux = flux,
      out = out,
      wave_lum = wave_lum,
      lum_unatten = lum_unatten,
      lum_atten = lum,
      lumtot_unatten = lumtot_unatten,
      lumtot_atten = lumtot_atten,
      lumtot_birth = lumtot_birth,
      lumtot_screen = lumtot_screen,
      agevec = agevec,
      SFR = SFR,
      masstot = masstot,
      massvec = massvec,
      M2L = masstot / lumtot_unatten,
      SFRburst = SFRburst,
      Zvec = Zvec,
      emission_input = emission_input
    )
  )
}

SMstarfunc = function(massfunc = massfunc_b5,
                      forcemass = FALSE,
                      agescale = 1,
                      burstage = c(0, 1e8),
                      youngage = c(1e8, 1e9),
                      midage = c(1e9, 5e9),
                      oldage = c(5e9, 9e9),
                      ancientage = c(9e9, 1.3e10),
                      stellpop = 'BC03lr',
                      speclib = NULL,
                      Z = 5,
                      z = 0.1,
                      H0 = 67.8,
                      OmegaM = 0.308,
                      OmegaL = 1 - OmegaM,
                      ref,
                      unimax = 13.8e9,
                      agemax = NULL,
                      intSFR = FALSE,
                      ...) {
  dots = list(...)
  massfunc_args = dots[names(dots) %in% names(formals(massfunc))]

  if (stellpop == 'BC03lr') {
    if (is.null(speclib)) {
      BC03lr = NULL
      data('BC03lr', envir = environment())
      speclib = BC03lr
    }
  }else if (stellpop == 'BC03hr') {
    if (is.null(speclib)) {
      BC03hr = NULL
      data('BC03hr', envir = environment())
      speclib = BC03hr
    }
  }else if (stellpop == 'EMILES') {
    if (is.null(speclib)) {
      EMILES = NULL
      data('EMILES', envir = environment())
      speclib = EMILES
    }
  }else if (stellpop == 'BPASS') {
    if (is.null(speclib)) {
      BPASS = NULL
      data('BPASS', envir = environment())
      speclib = BPASS
    }
  }

  if (any(speclib$Age <= 1e7)) {
    birthcloud = max(which(speclib$Age <= 1e7))
  } else{
    birthcloud = 1
  }

  if (!is.function(Z)) {
    if (Z %% 1 != 0) {
      #Check if the Z is non integer, if so then convert to a function
      tempZfunc = function(age, Z, ...) {
        rep(Z, length(age))
      }
      formals(tempZfunc)$Z = Z
      Z = tempZfunc
    }
  }

  agevec = speclib$Age * agescale

  if (is.function(Z)) {
    dots = list(...)
    Z_args = dots[names(dots) %in% names(formals(Z))]
    Zvec = do.call('Z', c(
      list(agevec),
      list(massfunc = massfunc),
      Z_args,
      massfunc_args
    ))
    Zlist = interp_param(Zvec, speclib$Z, log = TRUE)
    Zwmat = matrix(0, length(speclib$Age), length(speclib$Z))
    Zwmat[cbind(1:length(speclib$Age), Zlist[, 'ID_hi'])] = Zlist[, 'wt_hi']
    Zwmat[cbind(1:length(speclib$Age), Zlist[, 'ID_lo'])] = Zlist[, 'wt_lo']
    Zuse = which(colSums(Zwmat) > 0)
    Zdoweight = TRUE
  } else{
    Zuse = Z
    Zvec = rep(speclib$Z[Zuse], length(speclib$Age))
    Zdoweight = FALSE
  }

  if (intSFR) {
    massvec = rep(0, length(agevec))
    for (i in 1:length(agevec)) {
      tempint = try(do.call('integrate', c(
        list(
          f = massfunc,
          lower = speclib$AgeBins[i] * agescale,
          upper = speclib$AgeBins[i + 1] * agescale
        ),
        massfunc_args
      ))$value,
      silent = TRUE)
      if (inherits(tempint, "try-error")) {
        massvec[i] = 0
      } else{
        massvec[i] = tempint
      }
    }
    if (sum(massvec, na.rm = TRUE) == 0) {
      #pracma integral seems to work for very bursty star formation where integrate fails
      massvec = rep(0, length(agevec))
      for (i in 1:length(agevec)) {
        tempint = try(do.call('integral', c(
          list(
            f = massfunc,
            xmin = speclib$AgeBins[i] * agescale,
            xmax = speclib$AgeBins[i + 1] * agescale
          ),
          massfunc_args
        )),
        silent = TRUE)
        if (inherits(tempint, "try-error")) {
          massvec[i] = 0
        } else{
          massvec[i] = tempint
        }
      }
    }
  } else{
    massvec = do.call('massfunc', c(list(agevec), massfunc_args)) * speclib$AgeWeights
  }

  if (unimax != FALSE) {
    if (is.null(agemax)) {
      agemax = unimax - cosdistTravelTime(
        z = z,
        H0 = H0,
        OmegaM = OmegaM,
        OmegaL = OmegaL,
        ref = ref
      ) * 1e9
    }
  }

  if (!is.null(agemax)) {
    massvec[speclib$Age > agemax] = 0
  }

  if (forcemass == FALSE) {
    forcescale = 1
  } else{
    masstot = sum(massvec)
    forcescale = forcemass / masstot
    massvec = massvec * forcescale
  }

  totstar = rep(0, length(massvec))

  for (Zid in Zuse) {
    if (Zdoweight) {
      totstar = totstar + speclib$Zevo[[Zid]][, 'SMstar'] * massvec * Zwmat[, Zid]
    } else{
      totstar = speclib$Zevo[[Zid]][, 'SMstar'] * massvec
    }
  }

  burstageloc = c(min(which(speclib$Age - burstage[1] > 0)),
                  max(which(speclib$Age - burstage[2] < 0)))
  youngageloc = c(min(which(speclib$Age - youngage[1] > 0)),
                  max(which(speclib$Age - youngage[2] < 0)))
  midageloc = c(min(which(speclib$Age - midage[1] > 0)),
                max(which(speclib$Age - midage[2] < 0)))
  oldageloc = c(min(which(speclib$Age - oldage[1] > 0)),
                max(which(speclib$Age - oldage[2] < 0)))
  ancientageloc = c(min(which(speclib$Age - ancientage[1] > 0)),
                    max(which(speclib$Age - ancientage[2] < 0)))
  #burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  #youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  #midageloc=c(which.min(abs(speclib$Age-midage[1])),which.min(abs(speclib$Age-midage[2])))
  #oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  #ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))

  burstrescale = (burstage[2] - burstage[1]) / sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  youngrescale = (youngage[2] - youngage[1]) / sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  midrescale = (midage[2] - midage[1]) / sum(speclib$AgeWeights[midageloc[1]:midageloc[2]])
  oldrescale = (oldage[2] - oldage[1]) / sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientrescale = (ancientage[2] - ancientage[1]) / sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])

  burstform = sum(massvec[burstageloc[1]:burstageloc[2]]) * burstrescale
  burststar = sum(totstar[burstageloc[1]:burstageloc[2]]) * burstrescale
  youngform = sum(massvec[youngageloc[1]:youngageloc[2]]) * youngrescale
  youngstar = sum(totstar[youngageloc[1]:youngageloc[2]]) * youngrescale
  midform = sum(massvec[midageloc[1]:midageloc[2]]) * midrescale
  midstar = sum(totstar[midageloc[1]:midageloc[2]]) * midrescale
  oldform = sum(massvec[oldageloc[1]:oldageloc[2]]) * oldrescale
  oldstar = sum(totstar[oldageloc[1]:oldageloc[2]]) * oldrescale
  ancientform = sum(massvec[ancientageloc[1]:ancientageloc[2]]) * ancientrescale
  ancientstar = sum(totstar[ancientageloc[1]:ancientageloc[2]]) * ancientrescale


  return(
    c(
      BurstSMform = burstform,
      YoungSMform = youngform,
      MidSMform = midform,
      OldSMform = oldform,
      AncientSMform = ancientform,
      BurstSMstar = burststar,
      YoungSMstar = youngstar,
      MidSMstar = midstar,
      OldSMstar = oldstar,
      AncientSMstar = ancientstar,
      TotSMform = sum(massvec, na.rm = TRUE),
      TotSMstar = sum(totstar, na.rm = TRUE)
    )
  )
}

SFHburst = function(burstmass = 1e8,
                    burstage = 0,
                    stellpop = 'BC03lr',
                    speclib = NULL,
                    tau_birth = 1.0,
                    tau_screen = 0.3,
                    pow_birth = -0.7,
                    pow_screen = -0.7,
                    filters = 'all',
                    Z = 0.02,
                    emission = FALSE,
                    veldisp = 50,
                    range = 5,
                    res = 0.5,
                    emission_scale = 'FUV',
                    escape_frac = 1 - emission,
                    Ly_limit = 911.75,
                    LKL10 = NULL,
                    disp_stars = FALSE,
                    LSF = NULL,
                    z = 0.1,
                    H0 = 67.8,
                    OmegaM = 0.308,
                    OmegaL = 1 - OmegaM,
                    ref,
                    outtype = 'mag',
                    sparse = 5,
                    unimax = 13.8e9,
                    agemax = NULL,
                    LumDist_Mpc = NULL,
                    Eb = 0,
                    L0 = 2175.8,
                    LFWHM = 470,
                    ...) {
  burstmass = .interval(burstmass, 0, Inf, reflect = FALSE)

  if (stellpop == 'BC03lr') {
    if (is.null(speclib)) {
      BC03lr = NULL
      data('BC03lr', envir = environment())
      speclib = BC03lr
    }
  }else if (stellpop == 'BC03hr') {
    if (is.null(speclib)) {
      BC03hr = NULL
      data('BC03hr', envir = environment())
      speclib = BC03hr
    }
  }else if (stellpop == 'EMILES') {
    if (is.null(speclib)) {
      EMILES = NULL
      data('EMILES', envir = environment())
      speclib = EMILES
    }
  }else if (stellpop == 'BPASS') {
    if (is.null(speclib)) {
      BPASS = NULL
      data('BPASS', envir = environment())
      speclib = BPASS
    }
  }

  if (any(speclib$Age <= 1e7)) {
    birthcloud = max(which(speclib$Age <= 1e7))
  } else{
    birthcloud = 1
  }

  wave_lum = speclib$Wave
  if (sparse > 1) {
    sparse = seq(1, dim(speclib$Zspec[[1]])[2], by = sparse)
    wave_lum = wave_lum[sparse]
  } else{
    sparse = TRUE
  }

  if (!is.null(filters)) {
    if (filters[1] == 'all' | filters[1] == 'GAMA') {
      filters = c(
        'FUV_GALEX',
        'NUV_GALEX',
        'u_VST',
        'g_VST',
        'r_VST',
        'i_VST',
        'Z_VISTA',
        'Y_VISTA',
        'J_VISTA',
        'H_VISTA',
        'K_VISTA',
        'W1_WISE' ,
        'W2_WISE',
        'W3_WISE',
        'W4_WISE',
        'P100_Herschel',
        'P160_Herschel',
        'S250_Herschel' ,
        'S350_Herschel',
        'S500_Herschel'
      )
    }
  }

  if (unimax != FALSE) {
    if (is.null(agemax)) {
      agemax = unimax - cosdistTravelTime(
        z = z,
        H0 = H0,
        OmegaM = OmegaM,
        OmegaL = OmegaL,
        ref = ref
      ) * 1e9
    }
  }

  if (!is.null(agemax)) {
    if (burstage > agemax) {
      burstmass = 0
    }
  }

  metal_interp = interp_quick(Z, speclib$Z, log = TRUE)
  age_interp = interp_quick(burstage, speclib$Age)
  if (metal_interp['wt_lo'] == 1) {
    lum_unatten = speclib$Zspec[[metal_interp['ID_lo']]][age_interp['ID_lo'], sparse] *
      age_interp['wt_lo'] +
      speclib$Zspec[[metal_interp['ID_lo']]][age_interp['ID_hi'], sparse] *
      age_interp['wt_hi']
  } else{
    lum_unatten = speclib$Zspec[[metal_interp['ID_lo']]][age_interp['ID_lo'], sparse] *
      age_interp['wt_lo'] * metal_interp['wt_lo'] +
      speclib$Zspec[[metal_interp['ID_lo']]][age_interp['ID_hi'], sparse] *
      age_interp['wt_hi'] * metal_interp['wt_lo'] +
      speclib$Zspec[[metal_interp['ID_hi']]][age_interp['ID_lo'], sparse] *
      age_interp['wt_lo'] * metal_interp['wt_hi'] +
      speclib$Zspec[[metal_interp['ID_hi']]][age_interp['ID_hi'], sparse] *
      age_interp['wt_hi'] * metal_interp['wt_hi']
  }

  lum_unatten = lum_unatten * burstmass
  lum = lum_unatten
  lumtot_unatten = sum(.qdiff(wave_lum) * lum)

  if (tau_birth != 0 & burstage <= 1e7) {
    lum = lum * CF_birth(wave_lum, tau = tau_birth, pow = pow_birth)
    lumtot_birth = lumtot_unatten - sum(.qdiff(wave_lum) * lum)
  } else{
    lumtot_birth = 0
  }

  if(disp_stars){
    # Here we optionally disperse the spectrum along our LoS.

    grid = seq(-range, range, by=res)
    weights = dnorm(grid)
    wave_lum_log = log10(wave_lum)

    if(!is.null(LSF)){
      if(is.function(LSF)){
        vel_LSF = LSF(wave_lum*(1 + z)) #to get LSF dispersion in km/s into z frame
      }else if(is.matrix(LSF) | is.data.frame(LSF)){
        vel_LSF = approx(x=log10(LSF[,1]), y=LSF[,2], xout=log10(wave_lum*(1 + z)), rule=2)$y
      }else{
        vel_LSF = rep(LSF, length(wave_lum))
      }
    }else{
      vel_LSF = rep(0, length(wave_lum))
    }

    z_disp = sqrt(veldisp^2 + vel_LSF^2)/(.c_to_mps/1000) #this will be a vector of length wave_lum

    lum_conv = numeric(length(lum))

    for(i in seq_along(grid)){
      z_seq = grid[i]*z_disp
      new_wave = log10(wave_lum*(1 + z_seq))
      new_lum = log10(lum/(1 + z_seq))
      new_lum = 10^approx(x=new_wave, y=new_lum, xout=wave_lum_log, rule=2, yleft=new_lum[1], yright=new_lum[length(new_lum)])$y
      new_lum = new_lum*weights[i]
      lum_conv = lum_conv + new_lum
    }

    lum = lum_conv*res
    rm(lum_conv)
  }

  if (emission & burstage < 1e7) {
    if (emission_scale == 'FUV') {
      if (length(Ly_limit) == 1 | all(escape_frac == escape_frac[1])) {
        sel = which(wave_lum < Ly_limit[1])
        All_lum = (1 - escape_frac) * sum(.qdiff(wave_lum[sel]) * lum_unatten[sel])
      } else{
        All_lum = 0
        for (i in 1:(length(Ly_limit) - 1)) {
          sel = which(wave_lum < Ly_limit[i] & wave_lum > Ly_limit[i + 1])
          All_lum = All_lum + (1 - escape_frac[i]) * (Ly_limit[i] - Ly_limit[i + 1]) * mean(lum_unatten[sel])
        }
        sel = which(wave_lum < Ly_limit[length(Ly_limit)])
        All_lum = All_lum + (1 - escape_frac[length(Ly_limit)]) * sum(.qdiff(wave_lum[sel]) * lum_unatten[sel])
      }

      emission_input = list(All_lum = All_lum,
                            veldisp = veldisp,
                            Z = Z)
      emissionadd_unatten = emissionLines(All_lum = All_lum,
                                          veldisp = veldisp,
                                          Z = Z,
                                          LSF = LSF,
                                          z_LSF = z,
                                          LKL10 = LKL10,
                                          res = res,
                                          range = range)
    } else if (emission_scale == 'SFR') {
      SFRburst_emission = (1 - escape_frac) * burstmass / 1e7
      emission_input = list(SFR = SFRburst_emission,
                            veldisp = veldisp,
                            Z = Z)
      emissionadd_unatten = emissionLines(SFR = SFRburst_emission,
                                          veldisp = veldisp,
                                          Z = Z,
                                          LSF = LSF,
                                          z_LSF = z,
                                          LKL10 = LKL10,
                                          res = res,
                                          range = range)
    } else{
      stop('emission_scale must be one of SFR or FUV!')
    }

    emissionadd_atten = emissionadd_unatten
    emissionadd_atten$lum = emissionadd_atten$lum * CF_birth(emissionadd_atten$wave, tau = tau_birth, pow = pow_birth)

    lumtot_emission_unatten = sum(.qdiff(emissionadd_unatten$wave) * emissionadd_unatten$lum)
    lumtot_emission_atten = sum(.qdiff(emissionadd_atten$wave) * emissionadd_atten$lum)

    lumtot_birth = lumtot_birth - lumtot_emission_unatten + lumtot_emission_atten

    if(lumtot_birth < 0){
      lumtot_birth = 0
    }

    lum = addspec(wave_lum,
                  lum,
                  emissionadd_atten$wave,
                  emissionadd_atten$lum)
    lum_unatten = approxfun(log10(wave_lum), lum_unatten)(log10(lum$wave))

    wave_lum = lum$wave
    lum = lum$flux
    #wave_lum=lum_unatten$wave
    #lum_unatten=lum_unatten$flux
  } else{
    emission_input = NULL
  }

  if (tau_screen != 0) {
    lum = lum * CF_screen(
      wave_lum,
      tau = tau_screen,
      pow = pow_screen,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM
    )
    lumtot_screen = (lumtot_unatten - lumtot_birth) - sum(.qdiff(wave_lum) * lum)
  } else{
    lumtot_screen = 0
  }

  lumtot_atten = lumtot_unatten - sum(.qdiff(wave_lum) * lum)

  masstot = burstmass

  if (burstage <= 1e8) {
    SFRburst = masstot / 1e8
  } else{
    SFRburst = 0
  }

  if (z < 0 | is.null(filters)) {
    return(
      list(
        wave_lum = wave_lum,
        lum_atten = lum,
        lum_unatten = lum_unatten,
        lumtot_unatten = lumtot_unatten,
        lumtot_atten = lumtot_atten,
        lumtot_birth = lumtot_birth,
        lumtot_screen = lumtot_screen,
        agevec = burstage,
        SFR = burstmass,
        masstot = burstmass,
        M2L = burstmass / lumtot_unatten,
        SFRburst = SFRburst,
        Zvec = Z,
        emission_input = emission_input
      )
    ) # returns the minimal luminosity and mass outputs
  }
  if (z > 0) {
    flux = Lum2Flux(
      wave = wave_lum,
      lum = lum,
      z = z,
      H0 = H0,
      OmegaM = OmegaM,
      OmegaL = OmegaL,
      ref = ref,
      LumDist_Mpc = LumDist_Mpc
    )
    if (!is.null(outtype)) {
      out = photom_flux(flux, outtype = outtype, filters = filters)
      if (is.list(filters)) {
        cenout = {}
        for (i in filters) {
          cenout = c(cenout, cenwavefunc(i))
        }
        if (all(!is.null(names(filters)))) {
          out = data.frame(
            filter = names(filters),
            cenwave = cenout,
            out = out
          )
        } else{
          out = data.frame(filter = NA,
                           cenwave = cenout,
                           out = out)
        }
      } else{
        cenwave = NULL
        data('cenwave', envir = environment())
        out = data.frame(cenwave[match(filters, cenwave$filter),], out =
                           out)
      }
    } else{
      out = NULL
    }
  } else{
    flux = cbind(wave = wave_lum, flux = lum * .lsol_to_absolute)
    if (!is.null(outtype)) {
      out = photom_flux(flux, outtype = outtype, filters = filters)
      if (is.list(filters)) {
        cenout = {}
        for (i in filters) {
          cenout = c(cenout, cenwavefunc(i))
        }
        if (all(!is.null(names(filters)))) {
          out = data.frame(
            filter = names(filters),
            cenwave = cenout,
            out = out
          )
        } else{
          out = data.frame(filter = NA,
                           cenwave = cenout,
                           out = out)
        }
      } else{
        cenwave = NULL
        data('cenwave', envir = environment())
        out = data.frame(cenwave[match(filters, cenwave$filter),], out =
                           out)
      }
    } else{
      out = NULL
    }
  }

  return(
    list(
      flux = flux,
      out = out,
      wave_lum = wave_lum,
      lum_unatten = lum_unatten,
      lum_atten = lum,
      lumtot_unatten = lumtot_unatten,
      lumtot_atten = lumtot_atten,
      lumtot_birth = lumtot_birth,
      lumtot_screen = lumtot_screen,
      agevec = burstage,
      SFR = burstmass,
      masstot = burstmass,
      M2L = burstmass / lumtot_unatten,
      SFRburst = SFRburst,
      Zvec = Z,
      emission_input = emission_input
    )
  )
}
