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
                       addradio = FALSE,
                       Te = 1e4,
                       ff_frac = 0.1,
                       ff_power = -0.1,
                       sy_power = -0.8,
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
  call = match.call()
  
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
  
  
  if (is.null(AGN) | AGNlum == 0) {
    Final = SED_Stars_Bdust_Sdust
    AGN = NULL
    dustlum_AGN = 0
    dustmass_AGN = 0
  } else{
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
      AGN = AGN$final
      Final = data.frame(wave = SED_Stars_Bdust_Sdust$wave, flux = SED_Stars_Bdust_Sdust$flux +
                           AGN$flux)
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
  
  if (addradio) {
    Final = radiocont(
      wave = Final$wave,
      flux = Final$flux,
      z = 0,
      Te = Te,
      ff_frac = ff_frac,
      ff_power = ff_power,
      sy_power = sy_power,
      wavesamp = seq(6, max(waveout), by=0.1),
      flux_in = 'wave',
      flux_out = 'wave'
    )
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
    }
    filtout = list()
    for (i in filters) {
      filtout = c(filtout, list(approxfun(getfilt(i))))
    }
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
    photom_out = NULL
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
    DustEmit = data.frame(wave = Dust_Screen$Wave, lum = SED_Bdust_Sdust)
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
      call = call,
      z = z,
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
  
  Monitor = {
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
        list(filtout = quote(Data$filtout)),
        list(returnall = TRUE),
        list(Dale_M2L_func = quote(Data$Dale_M2L_func)),
        Data$arglist
      )
    )
    Photom = SEDout$Photom
    
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
        list(returnall = FALSE),
        Data$arglist
      )
    )
  }
  
  cutsig = (Data$flux$flux - Photom) / Data$flux$fluxerr
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
    names(Photom) = Data$flux$filter
    Monitor = c(Monitor, flux = Photom)
  }
  
  if (Data$fit == 'ld' | Data$fit == 'la' | Data$fit == 'check') {
    if ('LP' %in% Data$mon.names) {
      Monitor = c(Monitor, LP = LP)
    }
    if (length(Monitor) == 0) {
      Monitor = NA
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
        lwd = 5,
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
        lwd = 5,
        ...
      )
    }
    lines(x$StarsUnAtten,
          col = 'blue',
          lty = 2,
          lwd = 2)
    lines(x$StarsAtten, col = 'darkgreen', lwd = 2)
    lines(x$DustEmit, col = 'brown', lwd = 2)
    lines(x$AGN, col = 'purple', lwd = 2)
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
      lwd = c(5, 2, 2, 2, 2)
    )
    
    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = 5,
        grid = grid
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
      magicaxis::magaxis(4, col.axis = 'red', axis.col = 'red')
      legend(
        'bottomright',
        legend = c('SFR', 'Z'),
        col = c('black', 'red'),
        lty = 1,
        lwd = c(5, 2)
      )
    } else{
      plot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = 5
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
        lwd = c(5, 2)
      )
    }
  } else if (type == 'flux') {
    if (ylim[1] == 'auto') {
      ylim = quantile(x$FinalFlux[, 2], c(0.1, 1))
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
        lwd = 5,
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
        lwd = 5,
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
    points(x$Data$flux[, 2:3], pch = 16, col = 'red')
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magerr(x$Data$flux[, 2],
                        x$Data$flux[, 3],
                        ylo = x$Data$flux[, 4],
                        col = 'red')
    }
    legend('topleft', legend = paste('LP =', round(x$LP, 3)))
    
    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$Data$flux[, 2],
        (x$Data$flux[, 3] - x$SEDout$Photom) / x$Data$flux[, 4],
        pch = 16,
        col = 'red',
        grid = grid,
        log = 'x',
        xlim = xlim,
        ylim = c(-4, 4),
        xlab = xlab,
        ylab = '(Data-Model)/Error'
      )
    } else{
      plot(
        x$Data$flux[, 2],
        x$Data$flux[, 3] - x$SEDout$Photom,
        pch = 16,
        col = 'red',
        log = 'x',
        xlim = xlim,
        ylim = c(-4, 4),
        xlab = xlab,
        ylab = '(Data-Model)/Error'
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
