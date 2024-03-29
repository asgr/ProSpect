magAB2Jansky=function(x){10^(-0.4*(x-8.9))}
Jansky2magAB=function(x){-2.5*log10(x)+8.9}

magAB2CGS=function(x){10^(-0.4*(x+48.6))}
CGS2magAB=function(x){-2.5*log10(x)-48.6}

Jansky2CGS=function(x){x*1e-23}
CGS2Jansky=function(x){x*1e23}

AbsoluteToWHz=function(x){4*pi*((10*.pc_to_m)^2)*.jansky_to_si*magAB2Jansky(x)}
WHzToAbsolute=function(x){Jansky2magAB(x/(4*pi*((10*.pc_to_m)^2)*.jansky_to_si))}

magABcalc=function(wave, flux, filter='r_VST'){
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  if(is.character(filter)[1]){
    filter=getfilt(filter[1])
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter)
  return(-2.5 * log10(totlumnu) - 48.6)
}

CGScalc=function(wave, flux, filter='r_VST'){
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  if(is.character(filter)[1]){
    filter=getfilt(filter[1])
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter)
  return(totlumnu)
}

Janskycalc=function(wave, flux, filter='r_VST'){
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  if(is.character(filter)[1]){
    filter=getfilt(filter[1])
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter)
  return(totlumnu*1e23)
}


Lum2FluxFactor=function(z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, LumDist_Mpc=NULL){
  # Assuming lum to be converted is in the BC03 Lsol / Angstrom format
  # Because AB system is explicitly erg/s/cm^2/Hz flux
  if(z > 0){
    if(is.null(LumDist_Mpc)){
      LumDist_Mpc=cosdistLumDist(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref)
    }
    Dl_cm = LumDist_Mpc * .mpc_to_cm
    factor = .lsol_to_erg/(4*pi*Dl_cm^2)/(1 + z)
  }else{
    factor = .lsol_to_absolute
  }
  return(factor)
}

Lum2Flux=function(wave, lum, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, LumDist_Mpc=NULL){
  #Assumed lux input is Lsol / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  if(z > 0){
    if(is.null(LumDist_Mpc)){
      LumDist_Mpc = cosdistLumDist(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref)
    }
    Dl_cm = LumDist_Mpc * .mpc_to_cm
    flux = lum*.lsol_to_erg/(4*pi*Dl_cm^2)/(1 + z)
    wave = wave*(1 + z)
  }else{
    flux = lum * .lsol_to_absolute
  }
  #output is erg/s/cm^2/Ang (not per Hz! Need to make this final conversion to get to AB mag, but this is the standard way of viewing spectra).
  return(data.frame(wave=wave, flux=flux))
}

Flux2Lum=function(wave, flux, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, LumDist_Mpc=NULL){
  #Assumed flux input is erg/s/cm^2/Ang (not per Hz!)
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  if(z > 0){
    if(is.null(LumDist_Mpc)){
      LumDist_Mpc = cosdistLumDist(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref)
    }
    Dl_cm = LumDist_Mpc * .mpc_to_cm
    lum = flux/(.lsol_to_erg/(4*pi*Dl_cm^2)/(1 + z))
    wave = wave/(1 + z)
  }else{
    lum = flux/.lsol_to_absolute
  }
  #output is Lsol / Angstrom format
  return(data.frame(wave=wave, lum=lum))
}

photom_flux=function(wave, flux, outtype='mag', filters='all'){
  
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  
  if(is.matrix(filters) | is.data.frame(filters)){
    filters=list(filters)
  }
  
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
  
  if(outtype=='mag' | outtype=='magAB'){
    photom={}
    for(i in filters){
      photom=c(photom, magABcalc(wave=wave, flux=flux, filter=i))
    }
  }else if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    photom={}
    for(i in filters){
      photom=c(photom, Janskycalc(wave=wave, flux=flux, filter=i))
    }
  }else if(outtype=='cgs' | outtype=='CGS'){
    photom={}
    for(i in filters){
      photom=c(photom, CGScalc(wave=wave, flux=flux, filter=i))
    }
  }
  
  return(photom)
}

photom_lum=function(wave, lum, outtype='mag', filters='all', z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, LumDist_Mpc=NULL){
  
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  
  if(is.matrix(filters) | is.data.frame(filters)){
    filters=list(filters)
  }

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
  
  flux = Lum2Flux(wave=wave, lum=lum, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref, LumDist_Mpc=LumDist_Mpc)
    
  if(outtype=='mag' | outtype=='magAB'){
    photom={}
    for(i in filters){
      photom=c(photom, magABcalc(flux, filter=i))
    }
  }else if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    photom={}
    for(i in filters){
      photom=c(photom, Janskycalc(flux, filter=i))
    }
  }else if(outtype=='cgs' | outtype=='CGS'){
    photom={}
    for(i in filters){
      photom=c(photom, CGScalc(flux, filter=i))
    }
  }
  
  return(photom)
}

addspec=function(wave1, flux1, wave2, flux2, extrap='constant', waveout=NULL){
  wave1 = log10(wave1)
  wave2 = log10(wave2)
  if(anyNA(flux1)){
    flux1[is.na(flux1)] = 0
  }
  if(anyNA(flux2)){
    flux2[is.na(flux2)] = 0
  }
  flux1 = log10(flux1)
  flux2 = log10(flux2)

  if(is.null(waveout)){
    waveout = sort(unique(c(wave1,wave2)))
  }
  
  if(extrap=='constant'){
    flux1 = 10^approx(x=wave1, y=flux1, xout=waveout, rule=2, yleft=flux1[1], yright=flux1[length(flux1)])$y
    flux2 = 10^approx(x=wave2, y=flux2, xout=waveout, rule=2, yleft=flux2[1], yright=flux2[length(flux2)])$y
  }else{
    flux1 = 10^approx(x=wave1, y=flux1, xout=waveout, rule=2, yleft=log10(extrap), yright=log10(extrap))$y
    flux2 = 10^approx(x=wave2, y=flux2, xout=waveout, rule=2, yleft=log10(extrap), yright=log10(extrap))$y
  }

  # flux1 = 10^spline(x=wave1, y=flux1, xout=waveout)$y This behaves badly as discontinuities, so safer to use approx in general
  # flux2 = 10^spline(x=wave2, y=flux2, xout=waveout)$y This behaves badly as discontinuities, so safer to use approx in general
  
  flux1[is.na(flux1)] = 0
  flux2[is.na(flux2)] = 0
  return(data.frame(wave=10^waveout, flux=flux1+flux2))
}

#Lower level functions:

bandpass=function(wave, flux, filter, flux_in='freq', flux_out='freq', detect_type='photon'){
  # flux must be flux_nu, i.e. erg/s / cm^2 / Hz, not erg/s / cm^2 / Ang!
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  if(!is.function(filter)){
    filter = approxfun(x = filter[, 1], y = abs(filter[, 2]))
  }
  response = filter(wave)
  response[is.na(response)] = 0
  
  if(flux_in=='freq' & flux_out=='wave'){
    flux = convert_freq2wave(flux, wave)
    flux_in = 'wave'
  }
  
  if(flux_in=='wave' & flux_out=='freq'){
    flux = convert_wave2freq(flux, wave)
    flux_out = 'freq'
  }
  
  if(flux_out=='freq'){
    freq=1/wave
    freq_diff=abs(.qdiff(freq))
    if(detect_type=='photon'){
      output = response * wave * flux * freq_diff/sum(response * wave * freq_diff, na.rm = TRUE)
    }else if(detect_type=='energy'){
      output = response * flux * freq_diff/sum(response * freq_diff, na.rm = TRUE)
    }
  }else if(flux_out=='wave'){
    wave_diff=abs(.qdiff(wave))
    if(detect_type=='photon'){
      output = response * wave * flux * wave_diff/sum(response * wave * wave_diff, na.rm = TRUE)
    }else if(detect_type=='energy'){
      output = response * flux * wave_diff/sum(response * wave_diff, na.rm = TRUE)
    }
  }else{
    stop('flux_out must be one of: freq / wave.')
  }
  return(sum(output, na.rm=TRUE))
}

cenwavefunc=function(filter){
  if(is.function(filter)){
    return(NA)
  }else{
    wave=filter[,1]
    response=filter[,2]
    Ptot=sum(response, na.rm=TRUE)
    return((1/Ptot)*sum(response*wave, na.rm=TRUE))
  }
}

pivwavefunc=function(filter){
  if(is.function(filter)){
    return(NA)
  }else{
    wave=filter[,1]
    response=filter[,2]
    Ptot=sum(response/wave, na.rm=TRUE)
    return(sqrt((1/Ptot)*sum(response*wave, na.rm=TRUE)))
  }
}

convert_wave2freq=function(flux_wave, wave, wavefac=1e-10, freqfac=1){
  return(flux_wave*((wavefac/freqfac)*wave^2)/.c_to_mps)
}

convert_freq2wave=function(flux_freq, wave, wavefac=1e-10, freqfac=1){
  return(flux_freq*.c_to_mps/((wavefac/freqfac)*wave^2))
}
