specReBin = function(wave, flux, wavegrid=NULL, bin=NULL, binfunc=median,
                     interp='approx', logbin=TRUE, ...){
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  
  if(is.null(wavegrid)){
      if(logbin){
        wave = log10(wave)
        if(is.null(bin)){
          bin = binfunc(diff(wave))
        }
        wavegrid = seq(min(wave, na.rm=TRUE), max(wave, na.rm=TRUE), by=bin)
      }else{
        if(is.null(bin)){
          bin = binfunc(diff(wave))
        }
        wavegrid = seq(min(wave, na.rm=TRUE), max(wave, na.rm=TRUE), by=bin)
        wave = log10(wave)
        wavegrid = log10(wavegrid)
      }
  }else{
    wave = log10(wave)
    wavegrid = log10(wavegrid)
  }

  flux = log10(flux)
  
  if(interp == 'approx'){
    tempflux = approx(x=wave, y=flux, xout=wavegrid, ...)$y
  }else if(interp == 'spline'){
    tempflux = spline(x=wave, y=flux, xout=wavegrid, ...)$y
  }
  
  wavegrid = 10^wavegrid
  
  return(data.frame(wave=wavegrid, flux=tempflux))
}
