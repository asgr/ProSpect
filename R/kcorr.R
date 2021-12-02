kcorr = function(wave, lum, z, DistMod=NULL, filters='all', H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM, prospect=NULL){
  if(!is.null(prospect)){
    if(! inherits(prospect, "ProSpectSED")){stop('prospect input must be class ProSpectSED')}
    if(missing(wave)){wave = prospect$FinalLum$wave}
    if(missing(lum)){lum = prospect$FinalLum$lum}
    if(missing(z)){z = prospect$z}
    if(missing(filters)){filters = prospect$filtout}
  }
  
  absmag = photom_lum(wave=wave, lum=lum, z=0, outtype = 'magAB', filters=filters)
  apmag = photom_lum(wave=wave, lum=lum, z=z, outtype = 'magAB', filters=filters)
  
  if(is.null(DistMod)){
    DistMod = cosdistDistMod(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  
  kcorr = (apmag - DistMod) - absmag
  
  return(list(z=z, DistMod=DistMod, AbsMag=absmag, ApMag=apmag, kcorr=kcorr))
}

kcorr_evo = function(wave, lum, z=10^seq(-1,1,by=0.1), DistMod=NULL, filters='all', H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM, prospect=NULL){
  if(!is.null(prospect)){
    if(! inherits(prospect, "ProSpectSED")){stop('prospect input must be class ProSpectSED')}
    if(missing(wave)){wave = prospect$FinalLum$wave}
    if(missing(lum)){lum = prospect$FinalLum$lum}
    if(missing(filters)){filters = prospect$filtout}
  }
  
  absmag = photom_lum(wave=wave, lum=lum, z=0, outtype = 'magAB', filters=filters)
  apmag = {}
  
  if(is.null(DistMod)){
    DistMod = cosdistDistMod(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  
  kcorr = {}
  
  for(i in 1:length(DistMod)){
    apmag = rbind(apmag, photom_lum(wave=wave, lum=lum, z=z[i], outtype = 'magAB', filters=filters))
    kcorr = rbind(kcorr, (apmag[i,] - DistMod[i]) - absmag)
  }
  
  if(!is.null(names(filters))){
    colnames(apmag) = names(filters)
    colnames(kcorr) = names(filters)
  }
  
  return(list(z=z, DistMod=DistMod, AbsMag=absmag, ApMag=apmag, kcorr=kcorr))
}

Dcorr = function(wave, lum_atten, lum_unatten, filters='all', prospect=NULL){
  if(!is.null(prospect)){
    if(! inherits(prospect, "ProSpectSED")){stop('prospect input must be class ProSpectSED')}
    if(missing(wave)){wave = prospect$StarsAtten$wave}
    if(missing(lum_atten)){lum_atten = prospect$StarsAtten$lum}
    if(missing(lum_unatten)){lum_unatten = prospect$StarsUnAtten$lum}
    if(missing(filters)){filters = prospect$filtout}
  }
  
  absmag_atten = photom_lum(wave=wave, lum=lum_atten, z=0, outtype = 'magAB', filters=filters)
  absmag_unatten = photom_lum(wave=wave, lum=lum_unatten, z=0, outtype = 'magAB', filters=filters)
  
  Dcorr = absmag_atten - absmag_unatten
  
  return(list(AbsMag_atten = absmag_atten, AbsMag_unatten=absmag_unatten, Dcorr=Dcorr))
}
