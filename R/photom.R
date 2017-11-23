magABcalc=function(wave, flux, filter='r_VST', wavefac=1e-10){
  c=299792458
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  filter=getfilt(filter)
  fluxnu=(wavefac*flux*wave^2)/c
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
  return(-2.5 * log10(totlumnu) - 48.6)
}

CGScalc=function(wave, flux, filter='r_VST', wavefac=1e-10){
  c=299792458
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  filter=getfilt(filter)
  fluxnu=(wavefac*flux*wave^2)/c
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
  return(totlumnu)
}

Janskycalc=function(wave, flux, filter='r_VST', wavefac=1e-10){
  c=299792458
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  filter=getfilt(filter)
  fluxnu=(wavefac*flux*wave^2)/c
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
  return(totlumnu*1e23)
}


Lum2FluxFactor=function(z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref='planck'){
  #Assuming lum to be converted is in the BC03 Lsun / Angstrom format
  Lsun2erg= 3.826e33
  Mpc_to_cm = 3.08568e+26 # Because AB system is explicitly erg/s/cm^2/Angstrom flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)*Mpc_to_cm
  factor=Lsun2erg/(4*pi*Dl_cm^2)/(1+z)
  return(factor)
}

Lum2Flux=function(wave, lum, z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref='planck'){
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  #lum needs to be in the BC03 Lsun / Angstrom format
  Lsun2erg= 3.826e33
  Mpc_to_cm = 3.08568e+24 # Because AB system is explicitly erg/s/cm^2/Hz flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)*Mpc_to_cm
  flux=lum*Lsun2erg/(4*pi*Dl_cm^2)/(1+z)
  wave=wave*(1+z)
  #output is erg/s/cm^2/Angstrom (not per Hz! Need to make this final conversion to get to AB mag, but this is the standard way of viewing spectra).
  return(cbind(wave, flux))
}

photom=function(wave, lum, filters='all', z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref='planck'){
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  observedspec=Lum2Flux(wave=wave, lum = lum, z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)
  data('cenwave')
  if(filters[1]=='all'){filters=cenwave$filter}

  mag={}
  for(i in filters){
    mag=c(mag, magABcalc(observedspec, filter=i))
  }
  return(cbind(cenwave[match(filters, cenwave$filter),], mag=mag))
}
