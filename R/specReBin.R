specReBin = function(wave, flux, wavegrid=NULL, bin=NULL, binfunc=median,
                     interp='approx', logbin=TRUE, rough=TRUE, ...){
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
      tempflux = approx(x=wave, y=flux, xout=wavegrid_log, yleft=NA, yright=NA, ...)$y
    }else if(interp == 'spline'){
      tempflux = spline(x=wave, y=flux, xout=wavegrid_log, ...)$y
    }

    tempflux = 10^tempflux
  }else{
    # The input to spec_rebin_cpp always needs to be linear, even when logbin=TRUE
    # It computes the appropriate log bin limits internally
    tempflux = .spec_rebin_cpp(wave_in=wave, flux_in=flux, wave_out=wavegrid, logbin_in=logbin, logbin_out=logbin)
    
    # old code
    # if(interp == 'approx'){
    #   tempfun = approxfun(x=wave, y=flux, yleft=0, yright=0, ...)
    # }else if(interp == 'spline'){
    #   tempfun = spline(x=wave, y=flux, ...)
    # }
    # 
    # tempflux = rep(NA, length(wavegrid))
    # diffgrid = diff(wavegrid)
    # lendiff = length(diffgrid)
    # diffgrid = c(diffgrid[1], diffgrid)
    # for(i in 1:(lendiff + 1)){
    #   tempflux[i] = integrate(tempfun, wavegrid[i] - diffgrid[i]/2, wavegrid[i] + diffgrid[i]/2)$value/diffgrid[i]
    # }
  }
  
  return(data.frame(wave=wavegrid, flux=tempflux))
}
