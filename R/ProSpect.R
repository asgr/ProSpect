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
                          AGNct = 60,
                          AGNal = 4,
                          AGNbe = -0.5,
                          AGNta = 1,
                          AGNrm = 60,
                          AGNan = 30,
                          Eb = 0,
                          L0 = 2175.8,
                          LFWHM = 470,
                          IGMabsorb = 0,
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
  
  tau_birth = .interval(tau_birth, 0, 10, reflect = FALSE)
  tau_screen = .interval(tau_screen, 0, 10, reflect = FALSE)
  tau_AGN = .interval(tau_AGN, 0, 10, reflect = FALSE)
  pow_birth = .interval(pow_birth, -2, 0, reflect = FALSE)
  pow_screen = .interval(pow_screen, -2, 0, reflect = FALSE)
  pow_AGN = .interval(pow_AGN, -2, 0, reflect = FALSE)
  alpha_SF_birth = .interval(alpha_SF_birth, 0.0625, 4, reflect = FALSE)
  alpha_SF_screen = .interval(alpha_SF_screen, 0.0625, 4, reflect = FALSE)
  alpha_SF_AGN = .interval(alpha_SF_AGN, 0.0625, 4, reflect = FALSE)
  
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
  
  if (is.null(AGN) | AGNlum == 0) {
    #Final = Final
    AGN = NULL
    dustlum_AGN = 0
    dustmass_AGN = 0
  } else{
    if(isFALSE(Dale)){
      stop('Dale cannot be FALSE when using an AGN model!')
    }
    
    if (inherits(AGN, 'Fritz')) {
      #Use new model
      AGN = AGNinterp(
        lum = AGNlum,
        ct = AGNct,
        al = AGNal,
        be = AGNbe,
        ta = AGNta,
        rm = AGNrm,
        an = AGNan,
        Fritz = AGN
      )
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
          wavesamp = seq(6, waveout_max, by=0.1),
          flux_in = 'wave',
          flux_out = 'wave',
          subtractonly = !addradio_AGN # whether to add AGN radio or just subtract Dale radio
        )
      
      AGN = AGN$final
      if (length(Final$flux) == length(AGN$flux)) {
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
        wavesamp = seq(6, waveout_max, by=0.1),
        flux_in = 'wave',
        flux_out = 'wave',
        subtractonly = !addradio_AGN # whether to add AGN radio or just subtract Dale radio
      )
      
      AGN = AGN$final
      if (length(SED_Stars_Bdust_Sdust$flux) == length(AGN$flux)) {
        Final = data.frame(wave = SED_Stars_Bdust_Sdust$wave, flux = SED_Stars_Bdust_Sdust$flux +
                             AGN$flux)
      } else{
        Final = addspec(
          wave1 = SED_Stars_Bdust_Sdust$wave,
          flux1 = SED_Stars_Bdust_Sdust$flux,
          wave2 = AGN$wave,
          flux2 = AGN$flux
        )
      }
      colnames(AGN)[2] = 'lum'
    }
  }
  colnames(Final)[2] = 'lum'
  
  if (IGMabsorb > 0) {
    sel = which(Final$wave < 1215.67)
    Final$lum[sel] = Final$lum[sel] * (1 - IGMabsorb)
    sel = which(Final$wave < 911.75)
    Final$lum[sel] = 0
  }
  
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
  
  ## Implement Photoz fitting mode 
  sink_file_for_annoying_pracma_messages = tempfile()
  sink(sink_file_for_annoying_pracma_messages)
  if("z" %in% Data$parm.names & Data$arglist$photoz){
    
    ztest = parm["z"]
    if(Data$logged[Data$parm.names == "z"]){
      ztest = 10^ztest
    }
    
    if (!requireNamespace("celestial", quietly = TRUE)) {
      stop("The celestial package is needed for this function to work. Please install it from GitHub/ASGR", call. = FALSE)
    }
    
    agemax_new = celestial::cosdistUniAgeAtz(z = ztest, ref = Data$arglist$ref)*1e9 ##need to be in years 
    magemax_new = agemax_new/1e9 ## need to be in Gyr
    Zagemax_new = agemax_new/1e9
    LumDist_Mpc_new = celestial::cosdistLumDist(z = ztest, ref = Data$arglist$ref)
    
    ## Now update the args in Data
    Data$arglist$agemax = unname(agemax_new)
    Data$arglist$magemax = unname(magemax_new)
    Data$arglist$Zagemax = unname(Zagemax_new)
    Data$arglist$LumDist_Mpc = unname(LumDist_Mpc_new)
  }
  sink()
  unlink(sink_file_for_annoying_pracma_messages, recursive = T)
  
  Monitor = {}
  
  if (returnall) {
    SEDout = do.call(
      'ProSpectSED',
      args = c(
        parmlist,
        list(SFH = quote(Data$SFH)),
        list(speclib = quote(Data$speclib)),
        list(Dale = quote(Data$Dale)),
        list(AGN = quote(Data$AGN)),
        list(filtout = quote(Data$filtout)),
        list(filters = NULL),
        list(returnall = TRUE),
        list(Dale_M2L_func = quote(Data$Dale_M2L_func)),
        Data$arglist
      )
    )
    if(is.null(Data$filtout)){
      #this means we are in spec-z mode
      Photom = specReBin(wave = SEDout$Photom[,'wave'],
                         flux = SEDout$Photom[,'flux'],
                         wavegrid = Data$flux[,'wave'],
                         logbin = ifelse(is.null(Data$logbin), TRUE, Data$logbin),
                         rough = ifelse(is.null(Data$rough), TRUE, Data$rough)
      )[,'flux']
      SEDout$Photom = data.frame(wave = Data$flux[,'wave'],
                                 flux = Photom)
    }else{
      Photom = SEDout$Photom
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
        list(filtout = quote(Data$filtout)),
        list(filters = NULL),
        list(returnall = FALSE),
        Data$arglist
      )
    )
    
    if(is.null(Data$filtout)){
      #this means we are in spec-z mode
      Photom = specReBin(wave = Photom[,'wave'],
                        flux = Photom[,'flux'],
                        wavegrid = Data$flux[,'wave'],
                        logbin = ifelse(is.null(Data$logbin), TRUE, Data$logbin),
                        rough = ifelse(is.null(Data$rough), TRUE, Data$rough)
                        )[,'flux']
    }
  }
  
  cutsig = (Data$flux[,'flux'] - Photom) / Data$flux[,'fluxerr']
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
    names(Photom) = Data$flux[,'filter']
    Monitor = c(Monitor, flux = Photom)
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

plot.ProSpectSED = function(x,
                            xlim = c(1e3, 1e7),
                            ylim = 'auto',
                            xlab = 'Wavelength (Ang)',
                            ylab = 'auto',
                            grid = TRUE,
                            type = 'lum',
                            lwd_main = 5,
                            lwd_comp = 5,
                            ...) {
  if (type == 'lum') {
    if (ylim[1] == 'auto') {
      ylim = c(quantile(x$FinalLum[, 2], 0.45),
               max(x$StarsUnAtten, na.rm = TRUE))
    }
    if (ylab[1] == 'auto') {
      ylab = 'Luminosity Density (Lsol/Ang)'
    }
    layout(rbind(1, 2), heights = c(0.7, 0.3))
    par(oma = c(3.1, 3.1, 1.1, 2.1))
    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$FinalLum,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        ...
      )
    } else{
      plot(
        x$FinalLum,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        ...
      )
    }
    lines(x$StarsUnAtten,
          col = 'blue',
          lty = 2,
          lwd = lwd_comp)
    lines(x$StarsAtten, col = 'darkgreen', lwd = lwd_comp)
    lines(x$DustEmit, col = 'brown', lwd = lwd_comp)
    lines(x$AGN, col = 'purple', lwd = lwd_comp)
    legend(
      'topright',
      legend = c(
        'Total Lum',
        'Star Un-Atten',
        'Stars Atten',
        'Dust Emit',
        'AGN'
      ),
      col = c('black', 'blue', 'darkgreen', 'brown', 'purple'),
      lty = c(1, 2, 1, 1, 1),
      lwd = c(lwd_main, lwd_comp, lwd_comp, lwd_comp, lwd_comp)
    )
    
    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        majorn = c(5,3)
      )
      par(usr = c(
        par()$usr[1:2],
        -max(x$Stars$Zvec, na.rm = TRUE) * 0.04,
        max(x$Stars$Zvec, na.rm = TRUE) * 1.04
      ))
      lines(x$Stars$agevec / 1e9,
            x$Stars$Zvec,
            col = 'red',
            lwd = 2)
      magicaxis::magaxis(4, col.axis = 'red', axis.col = 'red', majorn=3)
      legend(
        'bottomright',
        legend = c('SFR', 'Z'),
        col = c('black', 'red'),
        lty = 1,
        lwd = c(lwd_main, lwd_comp)
      )
    } else{
      plot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = lwd_main
      )
      par(usr = c(
        par()$usr[1:2],
        -max(x$Stars$Zvec, na.rm = TRUE) * 0.04,
        max(x$Stars$Zvec, na.rm = TRUE) * 1.04
      ))
      lines(x$Stars$agevec / 1e9,
            x$Stars$Zvec,
            col = 'red',
            lwd = 2)
      axis(4, col = 'red', col.axis = 'red')
      legend(
        'bottomright',
        legend = c('SFR', 'Z'),
        col = c('black', 'red'),
        lty = 1,
        lwd = c(lwd_main, lwd_comp)
      )
    }
  } else if (type == 'flux') {
    if (ylim[1] == 'auto') {
      ylim = quantile(x$FinalFlux[, 'flux'], c(0.1, 1))
    }
    if (ylab[1] == 'auto') {
      ylab = 'Flux Density (Jansky)'
    }
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$FinalFlux,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        ...
      )
    } else{
      plot(
        x$FinalFlux,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        ...
      )
    }
  } else{
    stop('type argument must be one of lum or flux!')
  }
}

plot.ProSpectSEDlike = function(x,
                                xlim = c(1e3, 1e7),
                                ylim = 'auto',
                                xlab = 'Wavelength (Ang)',
                                ylab = 'auto',
                                grid = TRUE,
                                type = 'flux',
                                ...) {
  if (type == 'flux') {
    layout(rbind(1, 2), heights = c(0.7, 0.3))
    par(oma = c(3.1, 3.1, 1.1, 1.1))
    par(mar = c(0, 0, 0, 0))
    plot(
      x$SEDout,
      xlim = xlim,
      ylim = ylim,
      xlab = '',
      ylab = ylab,
      grid = grid,
      type = 'flux',
      ...
    )
    
    if(is.null(x$Data$filtout)){
      data_mode = 'spec'
      photom = x$SEDout$Photom[,'flux']
      comp_type = 'l'
    }else{
      data_mode = 'photom'
      photom = x$SEDout$Photom
      comp_type = 'p'
    }
    
    input_names = colnames(x$Data$flux)
    if('pivwave' %in% input_names){
      wavename = 'pivwave'
    }else if('cenwave' %in% input_names){
      wavename = 'cenwave'
    }else if('wave' %in% input_names){
      wavename = 'wave'
    }else{
      stop('wavelength column name not recognised, must')
    }
    
    if(data_mode == 'photom'){
      points(x$Data$flux[, c(wavename, 'flux')], pch = 16, col = 'red')
    }
    
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      if(data_mode == 'photom'){
        magicaxis::magerr(x$Data$flux[, wavename],
                          x$Data$flux[, 'flux'],
                          ylo = x$Data$flux[, 'fluxerr'],
                          col = 'red')
      }else{
        magicaxis::magerr(x$Data$flux[, wavename],
                          x$Data$flux[, 'flux'],
                          ylo = x$Data$flux[, 'fluxerr'],
                          col =  hsv(alpha=0.5),
                          poly = TRUE,
                          border = NA)
      }
    }
    
    legend('topleft', legend = paste('LP =', round(x$LP, 3)))
    
    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$Data$flux[, wavename],
        (x$Data$flux[, 'flux'] - photom) / x$Data$flux[, 'fluxerr'],
        type = comp_type,
        pch = 16,
        col = 'red',
        grid = grid,
        log = 'x',
        xlim = xlim,
        ylim = c(-4, 4),
        xlab = xlab,
        ylab = '(Data - Model)/Error'
      )
    } else{
      plot(
        x$Data$flux[, wavename],
        x$Data$flux[, 'flux'] - photom,
        type = comp_type,
        pch = 16,
        col = 'red',
        log = 'x',
        xlim = xlim,
        ylim = c(-4, 4),
        xlab = xlab,
        ylab = '(Data - Model)/Error'
      )
    }
  } else if (type == 'lum') {
    plot(
      x$SEDout,
      xlim = xlim,
      ylim = ylim,
      xlab = '',
      ylab = ylab,
      grid = grid,
      type = 'lum',
      ...
    )
  }
}
