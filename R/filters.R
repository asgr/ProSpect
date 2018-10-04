getfilt=function(filter){
  out=NA
  if(filter=='FUV_GALEX'){filt_FUV_GALEX=NULL; data('filt_FUV_GALEX', envir = environment()); out=filt_FUV_GALEX}
  if(filter=='NUV_GALEX'){filt_NUV_GALEX=NULL; data('filt_NUV_GALEX', envir = environment()); out=filt_NUV_GALEX}
  if(filter=='u_SDSS'){filt_u_SDSS=NULL; data('filt_u_SDSS', envir = environment()); out=filt_u_SDSS}
  if(filter=='g_SDSS'){filt_g_SDSS=NULL; data('filt_g_SDSS', envir = environment()); out=filt_g_SDSS}
  if(filter=='r_SDSS'){filt_r_SDSS=NULL; data('filt_r_SDSS', envir = environment()); out=filt_r_SDSS}
  if(filter=='i_SDSS'){filt_i_SDSS=NULL; data('filt_i_SDSS', envir = environment()); out=filt_i_SDSS}
  if(filter=='z_SDSS'){filt_z_SDSS=NULL; data('filt_z_SDSS', envir = environment()); out=filt_z_SDSS}
  if(filter=='u_VST'){filt_u_VST=NULL; data('filt_u_VST', envir = environment()); out=filt_u_VST}
  if(filter=='g_VST'){filt_g_VST=NULL; data('filt_g_VST', envir = environment()); out=filt_g_VST}
  if(filter=='r_VST'){filt_r_VST=NULL; data('filt_r_VST', envir = environment()); out=filt_r_VST}
  if(filter=='i_VST'){filt_i_VST=NULL; data('filt_i_VST', envir = environment()); out=filt_i_VST}
  if(filter=='z_VST'){filt_z_VST=NULL; data('filt_z_VST', envir = environment()); out=filt_z_VST}
  if(filter=='Z_VISTA'){filt_Z_VISTA=NULL; data('filt_Z_VISTA', envir = environment()); out=filt_Z_VISTA}
  if(filter=='Y_VISTA'){filt_Y_VISTA=NULL; data('filt_Y_VISTA', envir = environment()); out=filt_Y_VISTA}
  if(filter=='J_VISTA'){filt_J_VISTA=NULL; data('filt_J_VISTA', envir = environment()); out=filt_J_VISTA}
  if(filter=='H_VISTA'){filt_H_VISTA=NULL; data('filt_H_VISTA', envir = environment()); out=filt_H_VISTA}
  if(filter=='K_VISTA'){filt_K_VISTA=NULL; data('filt_K_VISTA', envir = environment()); out=filt_K_VISTA}
  if(filter=='Ks_VISTA'){filt_Ks_VISTA=NULL; data('filt_Ks_VISTA', envir = environment()); out=filt_Ks_VISTA}
  if(filter=='W1_WISE' | filter=='W1'){filt_W1_WISE=NULL; data('filt_W1_WISE', envir = environment()); out=filt_W1_WISE}
  if(filter=='W2_WISE' | filter=='W2'){filt_W2_WISE=NULL; data('filt_W2_WISE', envir = environment()); out=filt_W2_WISE}
  if(filter=='W3_WISE' | filter=='W3'){filt_W3_WISE=NULL; data('filt_W3_WISE', envir = environment()); out=filt_W3_WISE}
  if(filter=='W4_WISE' | filter=='W4'){filt_W4_WISE=NULL; data('filt_W4_WISE', envir = environment()); out=filt_W4_WISE}
  if(filter=='I1_Spitzer' | filter=='I1'){filt_I1_Spitzer=NULL; data('filt_I1_Spitzer', envir = environment()); out=filt_I1_Spitzer}
  if(filter=='I2_Spitzer' | filter=='I2'){filt_I2_Spitzer=NULL; data('filt_I2_Spitzer', envir = environment()); out=filt_I2_Spitzer}
  if(filter=='I3_Spitzer' | filter=='I3'){filt_I3_Spitzer=NULL; data('filt_I3_Spitzer', envir = environment()); out=filt_I3_Spitzer}
  if(filter=='I4_Spitzer' | filter=='I4'){filt_I4_Spitzer=NULL; data('filt_I4_Spitzer', envir = environment()); out=filt_I4_Spitzer}
  if(filter=='M24_Spitzer' | filter=='M24'){filt_M24_Spitzer=NULL; data('filt_M24_Spitzer', envir = environment()); out=filt_M24_Spitzer}
  if(filter=='M70_Spitzer' | filter=='M70'){filt_M70_Spitzer=NULL; data('filt_M70_Spitzer', envir = environment()); out=filt_M70_Spitzer}
  if(filter=='M160_Spitzer' | filter=='M160'){filt_M160_Spitzer=NULL; data('filt_M160_Spitzer', envir = environment()); out=filt_M160_Spitzer}
  if(filter=='P70_Herschel' | filter=='P70'){filt_P70_Herschel=NULL; data('filt_P70_Herschel', envir = environment()); out=filt_P70_Herschel}
  if(filter=='P100_Herschel' | filter=='P100'){filt_P100_Herschel=NULL; data('filt_P100_Herschel', envir = environment()); out=filt_P100_Herschel}
  if(filter=='P160_Herschel' | filter=='P160'){filt_P160_Herschel=NULL; data('filt_P160_Herschel', envir = environment()); out=filt_P160_Herschel}
  if(filter=='S250_Herschel' | filter=='S250'){filt_S250_Herschel=NULL; data('filt_S250_Herschel', envir = environment()); out=filt_S250_Herschel}
  if(filter=='S350_Herschel' | filter=='S350'){filt_S350_Herschel=NULL; data('filt_S350_Herschel', envir = environment()); out=filt_S350_Herschel}
  if(filter=='S500_Herschel' | filter=='S500'){filt_S500_Herschel=NULL; data('filt_S500_Herschel', envir = environment()); out=filt_S500_Herschel}
  if(filter=='S450_JCMT' | filter=='S450'){filt_S450_JCMT=NULL; data('filt_S450_JCMT', envir = environment()); out=filt_S450_JCMT}
  if(filter=='S850_JCMT' | filter=='S850'){filt_S850_JCMT=NULL; data('filt_S850_JCMT', envir = environment()); out=filt_S850_JCMT}
  return(out)
}

bandpass=function(wave, flux, filter, lum = TRUE){
  # flux must be flux_nu, i.e. erg/s / cm^2 / Hz, not erg/s / cm^2 / Ang!
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux=wave[,2]
      wave=wave[,1]
    }
  }
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

cenwavefunc=function(filter){
  wave=filter[,1]
  response=filter[,2]
  Ptot=sum(response)
  return((1/Ptot)*sum(response*wave))
}

convert_wave2freq=function(flux_wave, wave, wavefac=1e-10, freqfac=1){
  return=(wavefac*flux_wave*wave^2)/.c_to_mps
}

convert_freq2wave=function(flux_freq, wave, wavefac=1e-10, freqfac=1){
  return=flux_freq*.c_to_mps/(wavefac*wave^2)
}
