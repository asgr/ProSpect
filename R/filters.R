getfilt=function(filter){
  out=NA
  if(filter=='FUV'){data('filt_FUV_GALEX');out=filt_FUV_GALEX}
  if(filter=='NUV'){data('filt_NUV_GALEX');out=filt_NUV_GALEX}
  if(filter=='u_SDSS'){data('filt_u_SDSS');out=filt_u_SDSS}
  if(filter=='g_SDSS'){data('filt_g_SDSS');out=filt_g_SDSS}
  if(filter=='r_SDSS'){data('filt_r_SDSS');out=filt_r_SDSS}
  if(filter=='i_SDSS'){data('filt_i_SDSS');out=filt_i_SDSS}
  if(filter=='z_SDSS'){data('filt_z_SDSS');out=filt_z_SDSS}
  if(filter=='u_VST'){data('filt_u_VST');out=filt_u_VST}
  if(filter=='g_VST'){data('filt_g_VST');out=filt_g_VST}
  if(filter=='r_VST'){data('filt_r_VST');out=filt_r_VST}
  if(filter=='i_VST'){data('filt_i_VST');out=filt_i_VST}
  if(filter=='z_VST'){data('filt_z_VST');out=filt_z_VST}
  if(filter=='Z_VISTA'){data('filt_Z_VISTA');out=filt_Z_VISTA}
  if(filter=='Y_VISTA'){data('filt_Y_VISTA');out=filt_Y_VISTA}
  if(filter=='J_VISTA'){data('filt_J_VISTA');out=filt_J_VISTA}
  if(filter=='H_VISTA'){data('filt_H_VISTA');out=filt_H_VISTA}
  if(filter=='K_VISTA'){data('filt_K_VISTA');out=filt_K_VISTA}
  if(filter=='Ks_VISTA'){data('filt_Ks_VISTA');out=filt_Ks_VISTA}
  if(filter=='W1'){data('filt_W1_WISE');out=filt_W1_WISE}
  if(filter=='W2'){data('filt_W2_WISE');out=filt_W2_WISE}
  if(filter=='W3'){data('filt_W3_WISE');out=filt_W3_WISE}
  if(filter=='W4'){data('filt_W4_WISE');out=filt_W4_WISE}
  if(filter==100 | filter=='P100'){data('filt_P100_Herschel');out=filt_P100_Herschel}
  if(filter==160 | filter=='P160'){data('filt_P160_Herschel');out=filt_P160_Herschel}
  if(filter==250 | filter=='S250'){data('filt_S250_Herschel');out=filt_S250_Herschel}
  if(filter==350 | filter=='S350'){data('filt_S350_Herschel');out=filt_S350_Herschel}
  if(filter==450 | filter=='S450' | filter==500 | filter=='S500'){data('filt_S450_Herschel');out=filt_S450_Herschel}
  return(out)
}

bandpass=function(wave, flux, filter, lum = T){
  tempfunc = approxfun(x = filter[, 1], y = abs(filter[, 2]))
  tempremap = tempfunc(wave)
  tempremap[is.na(tempremap)] = 0
  if (lum) {
    return(sum(tempremap * wave * flux, na.rm = TRUE)/sum(tempremap * wave, na.rm = TRUE))
  }
  else {
    return(tempremap * wave * flux/sum(tempremap * wave, na.rm = TRUE))
  }
}

cenwavefunc=function(filt){
  wave=filt[,1]
  flux=filt[,2]
  Ptot=sum(flux)
  return((1/Ptot)*sum(flux*wave))
}
