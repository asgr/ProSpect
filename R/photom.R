magAB2Jansky=function(x){10^(-0.4*(x-8.9))}
Jansky2magAB=function(x){-2.5*log10(x)+8.9}

magAB2CGS=function(x){10^(-0.4*(x+48.6))}
CGS2magAB=function(x){-2.5*log10(x)-48.6}

Jansky2CGS=function(x){x*1e-23}
CGS2Jansky=function(x){x*1e23}

magABcalc=function(wave, flux, filter='r_VST'){
  #Data should be in erg/s / cm^2 / Angstrom
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
  filter=getfilt(filter)
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
  filter=getfilt(filter)
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
  filter=getfilt(filter)
  fluxnu=convert_wave2freq(flux, wave)
  totlumnu = bandpass(flux = fluxnu, wave = wave, filter = filter, lum = T)
  return(totlumnu*1e23)
}


Lum2FluxFactor=function(z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref='planck'){
  #Assuming lum to be converted is in the BC03 Lsun / Angstrom format
  # Because AB system is explicitly erg/s/cm^2/Angstrom flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)*.Mpc_to_cm
  factor=.Lsun2erg/(4*pi*Dl_cm^2)/(1+z)
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
  # Because AB system is explicitly erg/s/cm^2/Hz flux
  Dl_cm=cosdistLumDist(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)*.Mpc_to_cm
  flux=lum*.Lsun2erg/(4*pi*Dl_cm^2)/(1+z)
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
  flux=cbind(wave,flux)

  cenwave=NULL
  data('cenwave', envir = environment())
  
  if(filters[1]=='all'){filters=cenwave$filter}
  
  if(outtype=='mag' | outtype=='magAB'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, magABcalc(flux, filter=i))
    }
  }
  
  if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, Janskycalc(flux, filter=i))
    }
  }
  
  if(outtype=='cgs' | outtype=='CGS'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, CGScalc(flux, filter=i))
    }
  }
  
  return(cbind(cenwave[match(filters, cenwave$filter),], out=outdata))
}

photom_lum=function(wave, lum, outtype='mag', filters='all', z = 0.1, H0 = 100, OmegaM = 0.3, OmegaL = 1 - OmegaM - OmegaR, OmegaR = 0, w0 = -1, wprime = 0, ref='planck'){
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      lum=wave[,2]
      wave=wave[,1]
    }
  }
  flux=Lum2Flux(wave=wave, lum = lum, z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, OmegaR = OmegaR, w0 = w0, wprime = wprime, ref = ref)

  cenwave=NULL
  data('cenwave', envir = environment())
  
  if(filters[1]=='all'){filters=cenwave$filter}
  
  if(outtype=='mag' | outtype=='magAB'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, magABcalc(flux, filter=i))
    }
  }
  
  if(outtype=='jansky' | outtype=='Jansky' | outtype=='Jy'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, Janskycalc(flux, filter=i))
    }
  }
  
  if(outtype=='cgs' | outtype=='CGS'){
    outdata={}
    for(i in filters){
      outdata=c(outdata, CGScalc(flux, filter=i))
    }
  }
  
  return(cbind(cenwave[match(filters, cenwave$filter),], out=outdata))
}

addspec=function(wave1, flux1, wave2, flux2){
  wave=sort(c(wave1, wave2))
  return=cbind(wave=wave, flux=approxfun(wave1, flux1, rule=2)(wave)+approxfun(wave2, flux2, rule=2)(wave))
}
