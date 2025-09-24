ProSpectSED = function(SFH = SFHfunc,
                          z = 0.1,
                          tau_birth = 1,
                          tau_screen = 0.3,
                          tau_AGN = 1,
                          pow_birth = -0.7,
                          pow_screen = -0.7,
                          pow_AGN = -0.7,
                          alpha_SF_birth = 1,
                          alpha_SF_screen = 3,
                          alpha_SF_AGN = 0,
                          AGNlum = 0,
                          sparse = 5,
                          speclib = NULL,
                          Dale = NULL,
                          AGN = NULL,
                          filtout = NULL,
                          filters = 'all',
                          Dale_M2L_func = NULL,
                          returnall = TRUE,
                          H0 = 67.8,
                          OmegaM = 0.308,
                          OmegaL = 1 - OmegaM,
                          waveout = seq(2, 9.35, by = 0.01),
                          ref,
                          unimax = 13.8e9,
                          agemax = NULL,
                          LumDist_Mpc = NULL,
                          addradio_SF = FALSE,
                          addradio_AGN = FALSE,
                          Te_SF = 1e4,
                          ff_frac_SF = 0.1,
                          ff_power_SF = -0.1,
                          sy_power_SF = -0.8,
                          Te_AGN = 1e4,
                          ff_frac_AGN = 0.1,
                          ff_power_AGN = -0.1,
                          sy_power_AGN = -0.8,
                          AGNct = 40,
                          AGNrm = 60,
                          AGNan = 30,
                          AGNta = 1,
                          AGNal = 4,
                          AGNbe = -0.5,
                          AGNp = 1,
                          AGNq = 1,
                          Eb = 0,
                          L0 = 2175.8,
                          LFWHM = 470,
                          IGMabsorb = 0,
                          Inoue14_LAFcoef = NULL,
                          Inoue14_DLAcoef = NULL,
                          ...) {
  #call = match.call()

  if(!is.null(waveout)){
    waveout_max = max(waveout)
  }else{
    waveout_max = 9.35
  }

  if ('emission' %in% names(list(...))) {
    if (list(...)$emission & missing(waveout)) {
      waveout = NULL
    }
  }

  if(!is.null(SFH)){
    tau_birth = .interval(tau_birth, 0, 10, reflect = FALSE)
    tau_screen = .interval(tau_screen, 0, 10, reflect = FALSE)
    pow_birth = .interval(pow_birth, -2, 0, reflect = FALSE)
    pow_screen = .interval(pow_screen, -2, 0, reflect = FALSE)
    alpha_SF_birth = .interval(alpha_SF_birth, 0.0625, 4, reflect = FALSE)
    alpha_SF_screen = .interval(alpha_SF_screen, 0.0625, 4, reflect = FALSE)
  }

  if(!is.null(AGN)){
    tau_AGN = .interval(tau_AGN, 0, 10, reflect = FALSE)
    pow_AGN = .interval(pow_AGN, -2, 0, reflect = FALSE)
    alpha_SF_AGN = .interval(alpha_SF_AGN, 0.0625, 4, reflect = FALSE)
  }

  if(!is.null(SFH)){

    Stars = SFH(
      z = z,
      tau_birth = tau_birth,
      tau_screen = tau_screen,
      pow_birth = pow_birth,
      pow_screen = pow_screen,
      sparse = sparse,
      speclib = speclib,
      filters = NULL,
      unimax = unimax,
      agemax = agemax,
      Eb = Eb,
      L0 = L0,
      LFWHM = LFWHM,
      ...
    )

    if(!isFALSE(Dale)){
      Dust_Birth = Dale_interp(alpha_SF = alpha_SF_birth, Dale = Dale)
      Dust_Screen = Dale_interp(alpha_SF = alpha_SF_screen, Dale = Dale)

      SED_Bdust_Sdust = Dust_Birth$Aspec * Stars$lumtot_birth + Dust_Screen$Aspec *
        Stars$lumtot_screen
      SED_Stars_Bdust_Sdust = addspec(
        wave1 = Stars$wave_lum,
        flux1 = Stars$lum_atten,
        wave2 = Dust_Screen$Wave,
        flux2 = SED_Bdust_Sdust,
        extrap = 0,
        waveout = waveout
      )
    }else{
      Final = data.frame(wave = Stars$wave_lum, lum = Stars$lum_atten)
      Dust_Birth = NULL
      Dust_Screen = NULL
      SED_Bdust_Sdust = NULL
      SED_Stars_Bdust_Sdust = NULL
    }

    if (!is.null(Dale_M2L_func) & returnall) {
      dustlum_birth = Stars$lumtot_birth
      dustlum_screen = Stars$lumtot_screen
      dustmass_birth = Stars$lumtot_birth / Dale_M2L_func(alpha_SF_birth)
      dustmass_screen = Stars$lumtot_screen / Dale_M2L_func(alpha_SF_screen)
    } else{
      dustlum_birth = 0
      dustlum_screen = 0
      dustmass_birth = 0
      dustmass_screen = 0
    }

    ## adding radio contribution using just FIR generated from star-formation
    if (addradio_SF) {
      Final = radiocont(
        wave = SED_Stars_Bdust_Sdust$wave,
        flux = SED_Stars_Bdust_Sdust$flux,
        z = 0,
        Te = Te_SF,
        ff_frac = ff_frac_SF,
        ff_power = ff_power_SF,
        sy_power = sy_power_SF,
        wavesamp = seq(6, waveout_max, by=0.1),
        flux_in = 'wave',
        flux_out = 'wave'
      )
    }else if(!isFALSE(Dale)) {
      Final = SED_Stars_Bdust_Sdust
    }
  } else {
    Dust_Birth = NULL
    Dust_Screen = NULL
    SED_Bdust_Sdust = NULL
    SED_Stars_Bdust_Sdust = NULL
    dustlum_birth = NULL
    dustlum_screen = NULL
    dustmass_birth = NULL
    dustmass_screen = NULL
    Stars = NULL
    Final = NULL
  }

  if (is.null(AGN) | AGNlum == 0) {
    #Final = Final
    AGN = NULL
    dustlum_AGN = 0
    dustmass_AGN = 0
  } else {
    if(isFALSE(Dale)){
      stop('Dale cannot be FALSE when using an AGN model!')
    }

    if(inherits(AGN, 'Fritz') |
       inherits(AGN, 'SKIRTOR')) {
      #Use new model
      if (inherits(AGN, 'Fritz')) {
        AGN = Fritz_interp(
          lum = AGNlum,
          ct = AGNct,
          al = AGNal,
          be = AGNbe,
          ta = AGNta,
          rm = AGNrm,
          an = AGNan,
          Fritz = AGN
        )
      } else if (inherits(AGN, 'SKIRTOR')) {
        AGN = SKIRTOR_interp(
          lum = AGNlum,
          ta = AGNta,
          p = AGNp,
          q = AGNq,
          ct = AGNct,
          rm = AGNrm,
          an = AGNan,
          SKIRTOR = AGN
        )
      }

      dustlum_AGN = NA
      dustmass_AGN = NA
      AGN = atten_emit(
        wave = AGN$wave,
        flux = AGN$lum * .erg_to_lsol,
        tau = tau_screen,
        pow = pow_screen,
        alpha_SF = alpha_SF_screen,
        Dale = Dale,
        Dale_M2L_func = Dale_M2L_func,
        waveout = waveout,
        Eb = Eb,
        L0 = L0,
        LFWHM = LFWHM
      )
      if (!is.null(Dale_M2L_func) & returnall) {
        dustlum_screen = dustlum_screen + AGN$total_atten
        dustmass_screen = dustmass_screen + AGN$dustmass
      }

      ## subtracts off AGN contribution to the radio continuum unless you specifically request to add it back
      AGN$final = radiocont(
        wave = AGN$final$wave,
        flux = AGN$final$flux,
        z = 0,
        Te = Te_AGN,
        ff_frac = ff_frac_AGN,
        ff_power = ff_power_AGN,
        sy_power = sy_power_AGN,
        wavesamp = seq(6, waveout_max, by = 0.1),
        flux_in = 'wave',
        flux_out = 'wave',
        subtractonly = !addradio_AGN # whether to add AGN radio or just subtract Dale radio
      )

      AGN = AGN$final
      if (is.null(Final)) {
        Final = AGN
      } else if (length(Final$flux) == length(AGN$flux)) {
        Final = data.frame(wave = Final$wave, flux = Final$flux +
                             AGN$flux)
      } else {
        Final = addspec(
          wave1 = Final$wave,
          flux1 = Final$flux,
          wave2 = AGN$wave,
          flux2 = AGN$flux
        )
      }
      colnames(AGN)[2] = 'lum'
    } else{
      #Use old model
      #First we attenuate by the hot torus
      AGN = atten_emit(
        wave = AGN$Wave,
        flux = AGN$Aspec * AGNlum * .erg_to_lsol,
        tau = tau_AGN,
        pow = pow_AGN,
        alpha_SF = alpha_SF_AGN,
        Dale = Dale,
        Dale_M2L_func = Dale_M2L_func,
        waveout = waveout
      ) #no bump for the hot torus part

      if (!is.null(Dale_M2L_func) & returnall) {
        dustlum_AGN = AGN$total_atten
        dustmass_AGN = AGN$dustmass
      }
      #Second we re-attenuate the above by the screen (since it still has to pass out of the galaxy)
      AGN = atten_emit(
        wave = AGN$final$wave,
        flux = AGN$final$flux,
        tau = tau_screen,
        pow = pow_screen,
        alpha_SF = alpha_SF_screen,
        Dale = Dale,
        Dale_M2L_func = Dale_M2L_func,
        waveout = waveout,
        Eb = Eb,
        L0 = L0,
        LFWHM = LFWHM
      )

      if (!is.null(Dale_M2L_func) & returnall) {
        dustlum_screen = dustlum_screen + AGN$total_atten
        dustmass_screen = dustmass_screen + AGN$dustmass
      } else{
        dustlum_AGN = 0
        dustmass_AGN = 0
      }

      ## subtracts off AGN contribution to the radio continuum unless you specifically request to add it back
      AGN$final = radiocont(
        wave = AGN$final$wave,
        flux = AGN$final$flux,
        z = 0,
        Te = Te_AGN,
        ff_frac = ff_frac_AGN,
        ff_power = ff_power_AGN,
        sy_power = sy_power_AGN,
        wavesamp = seq(6, waveout_max, by = 0.1),
        flux_in = 'wave',
        flux_out = 'wave',
        subtractonly = !addradio_AGN # whether to add AGN radio or just subtract Dale radio
      )

      AGN = AGN$final
      if (is.null(Final)) {
        Final = AGN
      } else if (length(Final$flux) == length(AGN$flux)) {
        Final = data.frame(wave = Final$wave, flux = Final$flux +
                             AGN$flux)
      } else {
        Final = addspec(
          wave1 = Final$wave,
          flux1 = Final$flux,
          wave2 = AGN$wave,
          flux2 = AGN$flux
        )
      }
      colnames(AGN)[2] = 'lum'
    }
  }
  colnames(Final)[2] = 'lum'

  ## Previous ProSpect had IGM absorb done on the rest frame
  ## Now do to the observer frame later on
  # if (IGMabsorb > 0) {
  #   sel = which(Final$wave < 1215.67)
  #   Final$lum[sel] = Final$lum[sel] * (1 - IGMabsorb)
  #   sel = which(Final$wave < 911.75)
  #   Final$lum[sel] = 0
  # }

  if (is.null(filtout) & !is.null(filters)) {
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

    filtout = list()
    for (i in filters) {
      filtout = c(filtout, list(approxfun(getfilt(i))))
    }
    names(filtout) = filters
  }

  if (z > 0 & !is.null(filtout)) {
    Flux = Lum2Flux(
      wave = Final$wave,
      lum = Final$lum,
      z = z,
      H0 = H0,
      OmegaM = OmegaM,
      OmegaL = OmegaL,
      ref = ref,
      LumDist_Mpc = LumDist_Mpc
    )
    Flux$flux = convert_wave2freq(flux_wave = Flux$flux * .cgs_to_jansky,
                                  wave = Flux$wave)

    if (IGMabsorb > 0 & is.numeric(IGMabsorb)){
      sel = which(Flux$wave/(1+z) < 1215.67)
      Flux$flux[sel] = Flux$flux[sel] * (1 - IGMabsorb)
      sel = which(Flux$wave/(1+z) < 911.75)
      Flux$flux[sel] = 0
    }else if (IGMabsorb == "Inoue14"){
      Flux$flux = Flux$flux * Inoue14_IGM(Flux$wave, z, Inoue14_LAFcoef, Inoue14_DLAcoef)
    }

    photom_out = {}
    for (i in 1:length(filtout)) {
      photom_out = c(photom_out,
                     bandpass(
                       flux = Flux$flux,
                       wave = Flux$wave,
                       filter = filtout[[i]]
                     ))
    }
  } else if (z > 0 & is.null(filtout)) {
    Flux = Lum2Flux(
      wave = Final$wave,
      lum = Final$lum,
      z = z,
      H0 = H0,
      OmegaM = OmegaM,
      OmegaL = OmegaL,
      ref = ref,
      LumDist_Mpc = LumDist_Mpc
    )
    Flux$flux = convert_wave2freq(flux_wave = Flux$flux * .cgs_to_jansky,
                                  wave = Flux$wave)
    if (IGMabsorb > 0 & is.numeric(IGMabsorb)){
      sel = which(Flux$wave/(1+z) < 1215.67)
      Flux$flux[sel] = Flux$flux[sel] * (1 - IGMabsorb)
      sel = which(Flux$wave/(1+z) < 911.75)
      Flux$flux[sel] = 0
    }else if (IGMabsorb == "Inoue14"){
      Flux$flux = Flux$flux * Inoue14_IGM(Flux$wave, z, Inoue14_LAFcoef, Inoue14_DLAcoef)
    }

    photom_out = Flux
  } else if (z <= 0 & !is.null(filtout)) {
    Flux = cbind(wave = Final$wave,
                 flux = Final$lum * .lsol_to_absolute)
    photom_out = photom_flux(Flux, outtype = 'magAB', filters = filtout)
    #photom_out = photom_lum(Final, z=0, outtype = 'magAB', filters = filtout) same as above, but redundant
  } else{
    Flux = NULL
    photom_out = NULL
  }

  if (returnall) {
    StarsAtten = data.frame(wave = Stars$wave_lum, lum = Stars$lum_atten)
    StarsUnAtten = data.frame(wave = Stars$wave, lum = Stars$lum_unatten)

    if(!isFALSE(Dale)){
      DustEmit = data.frame(wave = Dust_Screen$Wave, lum = SED_Bdust_Sdust)
    }else{
      DustEmit = NULL
    }

    output = list(
      Photom = photom_out,
      FinalFlux = Flux,
      FinalLum = Final,
      StarsAtten = StarsAtten,
      StarsUnAtten = StarsUnAtten,
      DustEmit = DustEmit,
      AGN = AGN,
      Stars = Stars,
      dustmass = c(
        birth = dustmass_birth,
        screen = dustmass_screen,
        AGN = dustmass_AGN,
        total = sum(c(
          dustmass_birth, dustmass_screen, dustmass_AGN
        ), na.rm = TRUE)
      ),
      dustlum = c(
        birth = dustlum_birth,
        screen = dustlum_screen,
        AGN = dustlum_AGN,
        total = sum(c(
          dustlum_birth, dustlum_screen, dustlum_AGN
        ), na.rm = TRUE)
      ),
      #call = call,
      z = z,
      filters = filters,
      filtout = filtout,
      cosmo = list(H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
    )
    class(output) = 'ProSpectSED'
    return(output)
  } else{
    return(photom_out)
  }
}

ProSpectSEDlike = function(parm = c(8, 9, 10, 10, 0, -0.5, 0.2), Data) {
  if (is.null(Data$fit)) {
    Data$fit = 'optim'
  }
  Data$fit = tolower(Data$fit)
  if (is.null(Data$like)) {
    Data$like = 'st'
  }
  if (is.null(Data$mon.names)) {
    Data$mon.names = 'NULL'
  }
  if (is.null(Data$verbose)) {
    Data$verbose = FALSE
  }

  if (Data$fit == 'optim' | Data$fit == 'cma') {
    returnall = FALSE #fastest first!
  } else if (Data$fit == 'check') {
    returnall = TRUE
  } else if ((
    'masstot' %in% Data$mon.names |
    'SFRburst' %in% Data$mon.names |
    (length(grep(
      'dustmass', Data$mon.names
    )) > 0 |
    length(grep(
      'dustlum', Data$mon.names
    )) > 0) & (Data$fit == 'ld' | Data$fit == 'la')
  )) {
    returnall = TRUE
  } else{
    returnall = FALSE #just to be safe!
  }

  names(parm) = Data$parm.names

  if(is.null(Data$photom) & isTRUE(Data$mode == 'photom')){
    Data$photom = Data$flux
  }

  if(is.null(Data$spec) & (isTRUE(Data$mode == 'spec') | is.null(Data$filtout))){
    Data$spec = Data$flux
  }

  if(is.null(Data$photom) & isTRUE(Data$mode == 'both')){
    Data$photom = Data$flux
  }

  if((isTRUE(Data$mode == 'photom') | isTRUE(Data$mode == 'both')) & is.null(Data$filtout)){
    stop('Missing filout data!')
  }

  if((isTRUE(Data$mode == 'photom') | isTRUE(Data$mode == 'both')) & is.null(Data$photom)){
    stop('Missing photom data!')
  }

  if((isTRUE(Data$mode == 'spec') | isTRUE(Data$mode == 'both')) & is.null(Data$spec)){
    stop('Missing spec data!')
  }

  ## Implement Photoz fitting mode
  if(("photoz" %in% names(Data$arglist)) & ("z" %in% Data$parm.names)) {
    if (!requireNamespace("celestial", quietly = TRUE)) {
      stop("The celestial package is needed for this to work. Please install it from GitHub/ASGR", call. = FALSE)
    }
    if (!(c("ref") %in% names(Data$arglist))){
      if (!all(c("H0", "OmegaM", "OmegaL") %in% names(Data$arglist))){
        stop("Cosmology must be specified for photo-z")
      }
    }
    if(Data$arglist$photoz){

      z_genSF = Data$arglist$z_genSF ## What redshift should we start star formation?

      ztest = parm["z"]
      if(Data$logged[Data$parm.names == "z"]){
        ztest = 10^ztest
      }

      if(is.function(Data$arglist$IGMfunc)){
        Data$arglist$IGMabsorb = Data$arglist$IGMfunc(ztest)
      }else{
        if(is.null(Data$arglist$IGMfunc)){
          Data$arglist$IGMabsorb = 0 ## Default IGM absorption function
        }else if (Data$arglist$IGMfunc == "Songaila04"){
          Data$arglist$IGMabsorb = pnorm(ztest, mean = 3.8, sd = 1.2)
        }else if (Data$arglist$IGMfunc == "Inoue14"){
          Data$arglist$IGMabsorb = "Inoue14"
        }else if (is.numeric(Data$arglist$IGMfunc)){
          Data$arglist$IGMabsorb = Data$arglist$IGMfunc
        }else{
          stop("IGMfunc not valid type.")
        }
      }

      if(is.null(Data$arglist$ref)){
        agemax_new = (celestial::cosdistUniAgeAtz(ztest, H0 = Data$arglist$HO, OmegaM = Data$arglist$OmegaM, OmegaL = Data$arglist$OmegaL))*1e9 ##need to be in years
        if(!is.null(z_genSF)){
          agemax_new = (celestial::cosdistUniAgeAtz(ztest, H0 = Data$arglist$HO, OmegaM = Data$arglist$OmegaM, OmegaL = Data$arglist$OmegaL) - celestial::cosdistUniAgeAtz(z_genSF, H0 = Data$arglist$HO, OmegaM = Data$arglist$OmegaM, OmegaL = Data$arglist$OmegaL))*1e9 ##need to be in years
        }
        LumDist_Mpc_new = celestial::cosdistLumDist(z = ztest, H0 = Data$arglist$HO, OmegaM = Data$arglist$OmegaM, OmegaL = Data$arglist$OmegaL)
      }else{
        agemax_new = (celestial::cosdistUniAgeAtz(ztest, ref = Data$arglist$ref))*1e9 ##need to be in years
        if(!is.null(z_genSF)){
          agemax_new = (celestial::cosdistUniAgeAtz(ztest, ref = Data$arglist$ref) - celestial::cosdistUniAgeAtz(z_genSF, ref = Data$arglist$ref))*1e9 ##need to be in years
        }
        LumDist_Mpc_new = celestial::cosdistLumDist(z = ztest, ref = Data$arglist$ref)
      }

      magemax_new = agemax_new/1e9 ## need to be in Gyr
      Zagemax_new = agemax_new/1e9

      ## Now update the args in Data
      Data$arglist$agemax = unname(agemax_new)
      Data$arglist$magemax = unname(magemax_new)
      Data$arglist$Zagemax = unname(Zagemax_new)
      Data$arglist$LumDist_Mpc = unname(LumDist_Mpc_new)
    }
  }

  if (!is.null(Data$constraints)) {
    parm = Data$constraints(parm)
  }

  if (!is.null(Data$intervals)) {
    parm[parm < Data$intervals$lo] = Data$intervals$lo[parm < Data$intervals$lo]
    parm[parm > Data$intervals$hi] = Data$intervals$hi[parm > Data$intervals$hi]
  }

  if (!is.null(Data$logged)) {
    if (length(Data$logged) == 1) {
      if (Data$logged) {
        parmlist = 10 ^ parm
      } else{
        parmlist = parm
      }
    } else{
      parmlist = parm
      parmlist[Data$logged] = 10 ^ parm[Data$logged]
    }
  } else{
    parmlist = parm
  }

  if (Data$verbose) {
    print(parmlist)
  }

  Monitor = {}

  if(isTRUE(Data$mode == 'spec') | isTRUE(Data$mode == 'both')){
    #Because we will need to full spectrum for processing
    filtout = NULL
  }else{
    filtout = Data$filtout
  }

  if (returnall) {
    SEDout = do.call(
      'ProSpectSED',
      args = c(
        parmlist,
        list(SFH = quote(Data$SFH)),
        list(speclib = quote(Data$speclib)),
        list(Dale = quote(Data$Dale)),
        list(AGN = quote(Data$AGN)),
        list(filtout = quote(filtout)),
        list(filters = NULL),
        list(returnall = TRUE),
        list(Dale_M2L_func = quote(Data$Dale_M2L_func)),
        Data$arglist
      )
    )
    if(is.null(Data$filtout) | isTRUE(Data$mode == 'spec') | isTRUE(Data$mode == 'both')){
      #this means we are in spec-z mode
      if(isTRUE(Data$mode == 'both')){
        Phot_out = photom_flux(SEDout$Photom[,'wave'], convert_freq2wave(SEDout$Photom[,'flux'], SEDout$Photom[,'wave'])*.jansky_to_cgs,
                               outtype = 'Jansky', filters = Data$filtout)
      }

      if(!isTRUE(all.equal(SEDout$Photom[,'wave'], Data$spec[,'wave']))){ #if not already on the same grid we rebin
        Spec_out = specReBin(wave = SEDout$Photom[,'wave'],
                           flux = SEDout$Photom[,'flux'],
                           wavegrid = Data$spec[,'wave'],
                           logbin = ifelse(is.null(Data$logbin), TRUE, Data$logbin),
                           rough = ifelse(is.null(Data$rough), TRUE, Data$rough)
        )[,'flux']
        SEDout$Photom = data.frame(wave = Data$spec[,'wave'],
                                   flux = Spec_out)
      }else{
        Spec_out = SEDout$Photom[,'flux']
      }
    }else{
      Phot_out = SEDout$Photom
    }

    if (length(grep('dustmass', Data$mon.names)) > 0) {
      Monitor = c(dustmass = SEDout$dustmass)
    }
    if (length(grep('dustlum', Data$mon.names)) > 0) {
      Monitor = c(Monitor, dustlum = SEDout$dustlum)
    }
    if ('masstot' %in% Data$mon.names) {
      Monitor = c(Monitor, masstot = SEDout$Stars$masstot)
    }
    if ('SFRburst' %in% Data$mon.names) {
      Monitor = c(Monitor, SFRburst = SEDout$Stars$SFRburst)
    }
  } else{

    Photom = do.call(
      'ProSpectSED',
      args = c(
        parmlist,
        list(SFH = quote(Data$SFH)),
        list(speclib = quote(Data$speclib)),
        list(Dale = quote(Data$Dale)),
        list(AGN = quote(Data$AGN)),
        list(filtout = quote(filtout)),
        list(filters = NULL),
        list(returnall = FALSE),
        Data$arglist
      )
    )

    if(is.null(filtout)){
      if(is.null(Data$filtout) | isTRUE(Data$mode == 'spec') | isTRUE(Data$mode == 'both')){
        #this means we are in spec-z mode
        if(isTRUE(Data$mode == 'both')){
          Phot_out = photom_flux(Photom[,'wave'], convert_freq2wave(Photom[,'flux'], Photom[,'wave'])*.jansky_to_cgs,
                                 outtype = 'Jansky', filters = Data$filtout)
        }

        if(!isTRUE(all.equal(Photom[,'wave'], Data$spec[,'wave']))){ #if not already on the same grid we rebin
          Spec_out = specReBin(wave = Photom[,'wave'],
                            flux = Photom[,'flux'],
                            wavegrid = Data$spec[,'wave'],
                            logbin = ifelse(is.null(Data$logbin), TRUE, Data$logbin),
                            rough = ifelse(is.null(Data$rough), TRUE, Data$rough)
                            )[,'flux']
        }else{
          Spec_out = Photom[,'flux']
        }
      }
    }else{
      Phot_out = Photom
    }
  }

  if(!is.null(filtout)){
    #just in normal photometry mode
    cutsig = (Data$flux[,'flux'] - Phot_out) / Data$flux[,'fluxerr']
  }else{
    cutsig = (Data$spec[,'flux'] - Spec_out) / Data$spec[,'fluxerr']

    if(isTRUE(Data$mode == 'both')){ #special case
      cutsig = c(cutsig, (Data$photom[,'flux'] - Phot_out) / Data$photom[,'fluxerr'])
    }
  }

  if (Data$like == 'norm') {
    LL = sum(dnorm(x = cutsig, log = TRUE), na.rm = TRUE)
  } else if (Data$like == 'chisq') {
    LL = dchisq(sum(cutsig ^ 2),
                df = length(Data$filtout) - length(parm),
                log = TRUE)
  } else if (Data$like == 'st') {
    vardata = var(cutsig, na.rm = TRUE)
    dof = 2 * vardata / (vardata - 1)
    #dof=interval(dof,0,Inf)
    dof = max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    LL = as.numeric(sum(dt(
      cutsig, df = dof, log = TRUE
    ), na.rm = TRUE))
  } else{
    stop('Bad like option!')
  }

  if (is.null(Data$prior)) {
    LP = LL
  } else{
    LP = LL + as.numeric(Data$prior(parm))
  }
  if (Data$verbose) {
    print(LP)
  }

  if (length(grep('flux', Data$mon.names)) > 0) {
    names(Phot_out) = Data$flux[,'filter']
    Monitor = c(Monitor, flux = Phot_out)
  }

  if (Data$fit == 'ld' | Data$fit == 'la' | Data$fit == 'check') {
    if ('LP' %in% Data$mon.names) {
      Monitor = c(Monitor, LP = LP)
    }
    if (length(Monitor) == 0) {
      Monitor = 0L #this needs to be 0L is empty (not NA) Alas cannot return nothing!
    } else{
      Monitor = Monitor[match(Data$mon.names, names(Monitor))]
    }
  }

  # Various returns:

  if (Data$fit == 'optim' | Data$fit == 'cma') {
    return(-LP)
  } else if (Data$fit == 'ld' | Data$fit == 'la') {
    return(list(
      LP = LP,
      Dev = -2 * LL,
      Monitor = Monitor,
      yhat = 1,
      parm = parm
    ))
  } else if (Data$fit == 'check') {
    output = list(
      LP = LP,
      Dev = -2 * LL,
      Monitor = Monitor,
      yhat = 1,
      parm = parm,
      SEDout = SEDout,
      Data = Data
    )
    class(output) = 'ProSpectSEDlike'
    return(output)
  } else{
    return('Bad fit type!')
  }
}
