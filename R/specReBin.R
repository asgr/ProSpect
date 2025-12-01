specReBin = function(wave, flux, wavegrid=NULL, invar=NULL, bin=NULL, binfunc=median,
                     interp='approx', logbin=TRUE, rough=FALSE, ...){
  #Does not actually matter what units wave and flux are in as long as they are linear.

  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux = wave[,2]
      wave = wave[,1]
    }
  }

  if(is.null(wavegrid)){
    if(logbin){
      if(is.null(bin)){
        bin = binfunc(diff(log10(wave)))
      }
      wavegrid = 10^seq(log10(min(wave, na.rm=TRUE)), log10(max(wave, na.rm=TRUE)), by=bin)
    }else{
      if(is.null(bin)){
        bin = binfunc(diff(wave))
      }
      wavegrid = seq(min(wave, na.rm=TRUE), max(wave, na.rm=TRUE), by=bin)
    }
  }

  if(rough){
    #because it is always a bit better to work in log space when interpolating
    wave = log10(wave)
    wavegrid_log = log10(wavegrid)
    flux = log10(flux)

    if(interp == 'approx'){
      tempflux = approx(x=wave, y=flux, xout=wavegrid_log, ...)$y
    }else if(interp == 'spline'){
      tempflux = spline(x=wave, y=flux, xout=wavegrid_log, ...)$y
    }

    tempflux = 10^tempflux
  }else{
    # The input to spec_rebin_cpp always needs to be linear, even when logbin=TRUE
    # It computes the appropriate log bin limits internally
    tempflux = .spec_rebin_cpp(wave_in=wave, flux_in=flux, wave_out=wavegrid, invar_in=invar,
                               logbin_in=logbin, logbin_out=logbin)
  }

  return(data.frame(wave=wavegrid, flux=tempflux))
}

speclibReBin = function(speclib, wavegrid=NULL, bin=NULL, binfunc=median,
                        interp='approx', logbin=TRUE, rough=FALSE, cores=1L, check=FALSE, ...){

  doParallel::registerDoParallel(cores=cores)

  test = specReBin(speclib$Wave, speclib$Zspec[[1]][1,],
                   wavegrid=wavegrid, bin=bin, binfunc=binfunc,
                   interp=interp, logbin=logbin, rough=rough)

  wave_len = dim(test)[1]
  age_len = dim(speclib$Zspec[[1]])[1]

  i = NULL
  Zspec_rebin = foreach(i = 1:length(speclib$Zspec))%dopar%{
    tempmat = matrix(0, age_len, wave_len)
    for(j in 1:age_len){
      tempmat[j,] = specReBin(speclib$Wave, speclib$Zspec[[i]][j,],
                              wavegrid=test$wave, rough=rough)$flux
    }
    return(tempmat)
  }

  speclib$Zspec = Zspec_rebin
  speclib$Wave = test$wave

  if(check){
    speclib_check(speclib)
  }

  return(speclib)
}

speclibReGrid = function(speclib, logAge_steps = NULL, logZ_steps = NULL, Zsol=0.02, cores=8){
  #SSP_logAge_steps = log10(speclib$Age)
  #SSP_logZ_steps = log10(speclib$Z/Zsol)

  if(is.null(logAge_steps)){
    message('Using logAge grid from speclib!')
    logAge_steps = log10(speclib$Age)
    SSP_logAge_steps = logAge_steps
  }else{
    if(length(logAge_steps) > 1){
      message('Interpolating onto logAge grid different to speclib!')
    }
    SSP_logAge_steps = log10(speclib$Age)
    if(min(logAge_steps) < min(SSP_logAge_steps)){
      stop('logAge_steps minimum is less than minimum logAge present in the provided speclib!')
    }
    if(max(logAge_steps) > max(SSP_logAge_steps)){
      stop('logAge_steps maximum is more than maximum logAge present in the provided speclib!')
    }
  }

  SSP_logAge_steps[SSP_logAge_steps == -Inf] = -100
  logAge_steps[logAge_steps == -Inf] = -100

  if(is.null(logZ_steps)){
    message('Using logZ grid from speclib!')
    logZ_steps = log10(speclib$Z/Zsol)
    SSP_logZ_steps = logZ_steps
  }else{
    if(length(logZ_steps) > 1){
      message('Interpolating onto logZ grid different to speclib!')
    }
    SSP_logZ_steps = log10(speclib$Z/Zsol)
    if(min(logZ_steps) < min(SSP_logZ_steps)){
      stop('logZ_steps minimum is less than minimum logZ present in the provided speclib!')
    }
    if(max(logZ_steps) > max(SSP_logZ_steps)){
      stop('logZ_steps maximum is more than maximum logZ present in the provided speclib!')
    }
  }

  #
  logAge_step = NULL
  logZ_step = NULL

  if(length(logAge_steps) == 1 & length(logZ_steps) == 1){
    temp_Z = interp_quick(logZ_steps, SSP_logZ_steps)
    temp_Age = interp_quick(logAge_steps, SSP_logAge_steps)

    spec_Zlo_Alo = speclib$Zspec[[temp_Z['ID_lo']]][temp_Age['ID_lo'],]*temp_Z['wt_lo']*temp_Age['wt_lo']
    spec_Zlo_Ahi = speclib$Zspec[[temp_Z['ID_lo']]][temp_Age['ID_hi'],]*temp_Z['wt_lo']*temp_Age['wt_hi']
    spec_Zhi_Alo = speclib$Zspec[[temp_Z['ID_hi']]][temp_Age['ID_lo'],]*temp_Z['wt_hi']*temp_Age['wt_lo']
    spec_Zhi_Ahi = speclib$Zspec[[temp_Z['ID_hi']]][temp_Age['ID_hi'],]*temp_Z['wt_hi']*temp_Age['wt_hi']
    spec_out = spec_Zlo_Alo + spec_Zlo_Ahi + spec_Zhi_Alo + spec_Zhi_Ahi
    return(data.frame(wave = speclib$Wave, lum = spec_out))
  }else{
    cores = min(cores, length(logZ_steps), detectCores())
    registerDoParallel(cores=cores)
  }

  Zspec = foreach(logZ_step = logZ_steps)%dopar%{
    message('  ',logZ_step)
    temp_Z = interp_quick(logZ_step, SSP_logZ_steps)
    output = foreach(logAge_step = logAge_steps)%do%{
      #message('    ',logAge_step)
      temp_Age = interp_quick(logAge_step, SSP_logAge_steps)
      spec_Zlo_Alo = speclib$Zspec[[temp_Z['ID_lo']]][temp_Age['ID_lo'],]*temp_Z['wt_lo']*temp_Age['wt_lo']
      spec_Zlo_Ahi = speclib$Zspec[[temp_Z['ID_lo']]][temp_Age['ID_hi'],]*temp_Z['wt_lo']*temp_Age['wt_hi']
      spec_Zhi_Alo = speclib$Zspec[[temp_Z['ID_hi']]][temp_Age['ID_lo'],]*temp_Z['wt_hi']*temp_Age['wt_lo']
      spec_Zhi_Ahi = speclib$Zspec[[temp_Z['ID_hi']]][temp_Age['ID_hi'],]*temp_Z['wt_hi']*temp_Age['wt_hi']
      return(spec_Zlo_Alo + spec_Zlo_Ahi + spec_Zhi_Alo + spec_Zhi_Ahi)
    }
    output = do.call(rbind, output)
    return(as.matrix(output))
  }

  Zevo = foreach(logZ_step = logZ_steps)%dopar%{
    message('  ',logZ_step)
    temp_Z = interp_quick(logZ_step, SSP_logZ_steps)
    output = foreach(logAge_step = logAge_steps)%do%{
      #message('    ',logAge_step)
      temp_Age = interp_quick(logAge_step, SSP_logAge_steps)
      evo_Zlo_Alo = speclib$Zevo[[temp_Z['ID_lo']]][temp_Age['ID_lo'],]*temp_Z['wt_lo']*temp_Age['wt_lo']
      evo_Zlo_Ahi = speclib$Zevo[[temp_Z['ID_lo']]][temp_Age['ID_hi'],]*temp_Z['wt_lo']*temp_Age['wt_hi']
      evo_Zhi_Alo = speclib$Zevo[[temp_Z['ID_hi']]][temp_Age['ID_lo'],]*temp_Z['wt_hi']*temp_Age['wt_lo']
      evo_Zhi_Ahi = speclib$Zevo[[temp_Z['ID_hi']]][temp_Age['ID_hi'],]*temp_Z['wt_hi']*temp_Age['wt_hi']
      return(evo_Zlo_Alo + evo_Zlo_Ahi + evo_Zhi_Alo + evo_Zhi_Ahi)
    }
    output = as.data.frame(do.call(rbind, output))
    colnames(output) = colnames(speclib$Zevo[[1]])
    return(output)
  }

  Age_lims = .binlims(10^logAge_steps, log=T)
  AgeBins = c(Age_lims$lo[1], Age_lims$hi)
  AgeWeights = diff(AgeBins)

  logAge_steps[logAge_steps == -100] = -Inf

  SSP = list(
    Z = Zsol * 10^logZ_steps,
    Age = 10^logAge_steps,
    AgeBins = AgeBins,
    AgeWeights = AgeWeights,
    Wave = speclib$Wave,
    Labels = speclib$Labels,
    Zspec = Zspec,
    Zevo = Zevo
  )
}
