kcorr = function(wave, lum, z=NULL, DistMod=NULL, filters='all', H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM){
  absmag = photom_lum(wave=wave, lum=lum, z=0, outtype = 'magAB', filters=filters)
  apmag = photom_lum(wave=wave, lum=lum, z=z, outtype = 'magAB', filters=filters)
  
  if(is.null(DistMod)){
    DistMod = cosdistDistMod(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
  }
  
  kcorr = (apmag - DistMod) - absmag
  
  return(list(DistMod=DistMod, AbsMag=absmag, ApMag=apmag, kcorr=kcorr))
}

Dcorr = function(wave, lum_atten, lum_unatten, filters='all'){
  absmag_atten = photom_lum(wave=wave, lum=lum_atten, z=0, outtype = 'magAB', filters=filters)
  absmag_unatten = photom_lum(wave=wave, lum=lum_unatten, z=0, outtype = 'magAB', filters=filters)
  
  Dcorr = absmag_atten - absmag_unatten
  
  return(list(AbsMag_atten = absmag_atten, AbsMag_unatten=absmag_unatten, Dcorr=Dcorr))
}