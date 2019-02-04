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
  if(is.character(filter[1])){
    filter=getfilt(filter)
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
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
  if(is.character(filter[1])){
    filter=getfilt(filter)
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
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
  if(is.character(filter[1])){
    filter=getfilt(filter)
  }
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
  return(totlumnu*1e23)
}


Lum2FluxFactor=function(z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref){
  # Assuming lum to be converted is in the BC03 Lsun / Angstrom format
  # Because AB system is explicitly erg/s/cm^2/Angstrom flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref=ref)*.mpc_to_cm
  factor=.lsun_to_erg/(4*pi*Dl_cm^2)/(1+z)
  return(factor)
}

Lum2Flux=function(wave, lum, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref){
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  # Assuming lum to be converted is in the BC03 Lsun / Angstrom format
  # Because AB system is explicitly erg/s/cm^2/Hz flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref=ref)*.mpc_to_cm
  flux=lum*.lsun_to_erg/(4*pi*Dl_cm^2)/(1+z)
  wave=wave*(1+z)
  #output is erg/s/cm^2/Angstrom (not per Hz! Need to make this final conversion to get to AB mag, but this is the standard way of viewing spectra).
  return(cbind(wave=wave, flux=flux))
}

photom_flux=function(wave, flux, outtype='mag', filters='all'){
  
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  
  if(filters[1]=='all'){
    cenwave=NULL
    data('cenwave', envir = environment())
    filters=cenwave$filter
  }
  
  if(outtype=='mag' | outtype=='magAB'){
    photom={}
    for(i in filters){
      photom=c(photom, magABcalc(wave=wave, flux=flux, filter=i))
    }
  }
  
  if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    photom={}
    for(i in filters){
      photom=c(photom, Janskycalc(wave=wave, flux=flux, filter=i))
    }
  }
  
  if(outtype=='cgs' | outtype=='CGS'){
    photom={}
    for(i in filters){
      photom=c(photom, CGScalc(wave=wave, flux=flux, filter=i))
    }
  }
  
  return(photom)
}

photom_lum=function(wave, lum, outtype='mag', filters='all', z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref){
  
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  
  flux=Lum2Flux(wave=wave, lum=lum, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref)

  if(filters[1]=='all'){
    cenwave=NULL
    data('cenwave', envir = environment())
    filters=cenwave$filter
  }
  
  if(outtype=='mag' | outtype=='magAB'){
    photom={}
    for(i in filters){
      photom=c(photom, magABcalc(flux, filter=i))
    }
  }
  
  if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    photom={}
    for(i in filters){
      photom=c(photom, Janskycalc(flux, filter=i))
    }
  }
  
  if(outtype=='cgs' | outtype=='CGS'){
    photom={}
    for(i in filters){
      photom=c(photom, CGScalc(flux, filter=i))
    }
  }
  
  return(photom)
}

addspec=function(wave1, flux1, wave2, flux2){
  wave1=log10(wave1)
  wave2=log10(wave2)
  flux1=log10(flux1)
  flux2=log10(flux2)
  wave=sort(c(wave1, wave2))
  return=cbind(wave=10^wave, flux=10^approxfun(wave1, flux1, rule=2, yleft=flux1[1], yright=flux1[length(flux1)])(wave)+10^approxfun(wave2, flux2, rule=2, , yleft=flux2[1], yright=flux2[length(flux2)])(wave))
}
