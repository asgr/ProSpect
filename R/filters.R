.tophat = function(wavelo, wavehi){
  return(cbind(wave=c(wavelo,wavehi), response=c(1,1)))
}

getfilt=function(filter){
  out = NULL
  if(filter=='FUV_GALEX' | filter=='FUV'){filt_FUV_GALEX=NULL; data('filt_FUV_GALEX', envir = environment()); out = filt_FUV_GALEX}
  if(filter=='NUV_GALEX' | filter=='NUV'){filt_NUV_GALEX=NULL; data('filt_NUV_GALEX', envir = environment()); out = filt_NUV_GALEX}

  if(filter=='u_SDSS' | filter=='u'){filt_u_SDSS=NULL; data('filt_u_SDSS', envir = environment()); out = filt_u_SDSS}
  if(filter=='g_SDSS' | filter=='g'){filt_g_SDSS=NULL; data('filt_g_SDSS', envir = environment()); out = filt_g_SDSS}
  if(filter=='r_SDSS' | filter=='r'){filt_r_SDSS=NULL; data('filt_r_SDSS', envir = environment()); out = filt_r_SDSS}
  if(filter=='i_SDSS' | filter=='i'){filt_i_SDSS=NULL; data('filt_i_SDSS', envir = environment()); out = filt_i_SDSS}
  if(filter=='z_SDSS' | filter=='z'){filt_z_SDSS=NULL; data('filt_z_SDSS', envir = environment()); out = filt_z_SDSS}

  if(filter=='u_VST'){filt_u_VST=NULL; data('filt_u_VST', envir = environment()); out = filt_u_VST}
  if(filter=='g_VST'){filt_g_VST=NULL; data('filt_g_VST', envir = environment()); out = filt_g_VST}
  if(filter=='r_VST'){filt_r_VST=NULL; data('filt_r_VST', envir = environment()); out = filt_r_VST}
  if(filter=='i_VST'){filt_i_VST=NULL; data('filt_i_VST', envir = environment()); out = filt_i_VST}
  if(filter=='z_VST'){filt_z_VST=NULL; data('filt_z_VST', envir = environment()); out = filt_z_VST}

  if(filter=='u_LSST' | filter=='u_Rubin'){filt_u_LSST=NULL; data('filt_u_LSST', envir = environment()); out = filt_u_LSST}
  if(filter=='g_LSST' | filter=='r_Rubin'){filt_g_LSST=NULL; data('filt_g_LSST', envir = environment()); out = filt_g_LSST}
  if(filter=='r_LSST' | filter=='r_Rubin'){filt_r_LSST=NULL; data('filt_r_LSST', envir = environment()); out = filt_r_LSST}
  if(filter=='i_LSST' | filter=='i_Rubin'){filt_i_LSST=NULL; data('filt_i_LSST', envir = environment()); out = filt_i_LSST}
  if(filter=='z_LSST' | filter=='z_Rubin'){filt_z_LSST=NULL; data('filt_z_LSST', envir = environment()); out = filt_z_LSST}
  if(filter=='y_LSST' | filter=='Y_LSST' | filter=='y_Rubin' | filter=='Y_Rubin'){filt_y_LSST=NULL; data('filt_y_LSST', envir = environment()); out = filt_y_LSST}

  if(filter=='g_HSC'){filt_g_HSC=NULL; data('filt_g_HSC', envir = environment()); out = filt_g_HSC}
  if(filter=='r_HSC'){filt_r_HSC=NULL; data('filt_r_HSC', envir = environment()); out = filt_r_HSC}
  if(filter=='i_HSC'){filt_i_HSC=NULL; data('filt_i_HSC', envir = environment()); out = filt_i_HSC}
  if(filter=='z_HSC'){filt_z_HSC=NULL; data('filt_z_HSC', envir = environment()); out = filt_z_HSC}
  if(filter=='y_HSC' | filter=='Y_HSC'){filt_Y_HSC=NULL; data('filt_Y_HSC', envir = environment()); out = filt_Y_HSC}

  if(filter=='Z_VISTA' | filter=='Z'){filt_Z_VISTA=NULL; data('filt_Z_VISTA', envir = environment()); out = filt_Z_VISTA}
  if(filter=='Y_VISTA' | filter=='Y'){filt_Y_VISTA=NULL; data('filt_Y_VISTA', envir = environment()); out = filt_Y_VISTA}
  if(filter=='J_VISTA' | filter=='J'){filt_J_VISTA=NULL; data('filt_J_VISTA', envir = environment()); out = filt_J_VISTA}
  if(filter=='H_VISTA' | filter=='H'){filt_H_VISTA=NULL; data('filt_H_VISTA', envir = environment()); out = filt_H_VISTA}
  if(filter=='K_VISTA' | filter=='K'){filt_K_VISTA=NULL; data('filt_K_VISTA', envir = environment()); out = filt_K_VISTA}
  if(filter=='Ks_VISTA' | filter=='Ks'){filt_Ks_VISTA=NULL; data('filt_Ks_VISTA', envir = environment()); out = filt_Ks_VISTA}

  if(filter=='Z_UKIRT'){filt_Z_UKIRT=NULL; data('filt_Z_UKIRT', envir = environment()); out = filt_Z_UKIRT}
  if(filter=='Y_UKIRT'){filt_Y_UKIRT=NULL; data('filt_Y_UKIRT', envir = environment()); out = filt_Y_UKIRT}
  if(filter=='J_UKIRT'){filt_J_UKIRT=NULL; data('filt_J_UKIRT', envir = environment()); out = filt_J_UKIRT}
  if(filter=='H_UKIRT'){filt_H_UKIRT=NULL; data('filt_H_UKIRT', envir = environment()); out = filt_H_UKIRT}
  if(filter=='K_UKIRT'){filt_K_UKIRT=NULL; data('filt_K_UKIRT', envir = environment()); out = filt_K_UKIRT}

  if(filter=='F070W_JWST'){filt_F070W_JWST=NULL; data('filt_F070W_JWST', envir = environment()); out = filt_F070W_JWST}
  if(filter=='F090W_JWST'){filt_F090W_JWST=NULL; data('filt_F090W_JWST', envir = environment()); out = filt_F090W_JWST}
  if(filter=='F115W_JWST'){filt_F115W_JWST=NULL; data('filt_F115W_JWST', envir = environment()); out = filt_F115W_JWST}
  if(filter=='F140M_JWST'){filt_F140M_JWST=NULL; data('filt_F140M_JWST', envir = environment()); out = filt_F140M_JWST}
  if(filter=='F150W_JWST'){filt_F150W_JWST=NULL; data('filt_F150W_JWST', envir = environment()); out = filt_F150W_JWST}
  if(filter=='F162M_JWST'){filt_F162M_JWST=NULL; data('filt_F162M_JWST', envir = environment()); out = filt_F162M_JWST}
  if(filter=='F164N_JWST'){filt_F164N_JWST=NULL; data('filt_F164N_JWST', envir = environment()); out = filt_F164N_JWST}
  if(filter=='F150W2_JWST'){filt_F150W2_JWST=NULL; data('filt_F150W2_JWST', envir = environment()); out = filt_F150W2_JWST}
  if(filter=='F182M_JWST'){filt_F182M_JWST=NULL; data('filt_F182M_JWST', envir = environment()); out = filt_F182M_JWST}
  if(filter=='F187N_JWST'){filt_F187N_JWST=NULL; data('filt_F187N_JWST', envir = environment()); out = filt_F187N_JWST}
  if(filter=='F200W_JWST'){filt_F200W_JWST=NULL; data('filt_F200W_JWST', envir = environment()); out = filt_F200W_JWST}
  if(filter=='F210M_JWST'){filt_F210M_JWST=NULL; data('filt_F210M_JWST', envir = environment()); out = filt_F210M_JWST}
  if(filter=='F212N_JWST'){filt_F212N_JWST=NULL; data('filt_F212N_JWST', envir = environment()); out = filt_F212N_JWST}
  if(filter=='F250M_JWST'){filt_F250M_JWST=NULL; data('filt_F250M_JWST', envir = environment()); out = filt_F250M_JWST}
  if(filter=='F277W_JWST'){filt_F277W_JWST=NULL; data('filt_F277W_JWST', envir = environment()); out = filt_F277W_JWST}
  if(filter=='F300M_JWST'){filt_F300M_JWST=NULL; data('filt_F300M_JWST', envir = environment()); out = filt_F300M_JWST}
  if(filter=='F323N_JWST'){filt_F323N_JWST=NULL; data('filt_F323N_JWST', envir = environment()); out = filt_F323N_JWST}
  if(filter=='F322W2_JWST'){filt_F322W2_JWST=NULL; data('filt_F322W2_JWST', envir = environment()); out = filt_F322W2_JWST}
  if(filter=='F335M_JWST'){filt_F335M_JWST=NULL; data('filt_F335M_JWST', envir = environment()); out = filt_F335M_JWST}
  if(filter=='F356W_JWST'){filt_F356W_JWST=NULL; data('filt_F356W_JWST', envir = environment()); out = filt_F356W_JWST}
  if(filter=='F360M_JWST'){filt_F360M_JWST=NULL; data('filt_F360M_JWST', envir = environment()); out = filt_F360M_JWST}
  if(filter=='F405N_JWST'){filt_F405N_JWST=NULL; data('filt_F405N_JWST', envir = environment()); out = filt_F405N_JWST}
  if(filter=='F410M_JWST'){filt_F410M_JWST=NULL; data('filt_F410M_JWST', envir = environment()); out = filt_F410M_JWST}
  if(filter=='F430M_JWST'){filt_F430M_JWST=NULL; data('filt_F430M_JWST', envir = environment()); out = filt_F430M_JWST}
  if(filter=='F444W_JWST'){filt_F444W_JWST=NULL; data('filt_F444W_JWST', envir = environment()); out = filt_F444W_JWST}
  if(filter=='F460M_JWST'){filt_F460M_JWST=NULL; data('filt_F460M_JWST', envir = environment()); out = filt_F460M_JWST}
  if(filter=='F466N_JWST'){filt_F466N_JWST=NULL; data('filt_F466N_JWST', envir = environment()); out = filt_F466N_JWST}
  if(filter=='F470N_JWST'){filt_F470N_JWST=NULL; data('filt_F470N_JWST', envir = environment()); out = filt_F470N_JWST}
  if(filter=='F480M_JWST'){filt_F480M_JWST=NULL; data('filt_F480M_JWST', envir = environment()); out = filt_F480M_JWST}

  if(filter=='F560W_JWST'){filt_F560W_JWST=NULL; data('filt_F560W_JWST', envir = environment()); out = filt_F560W_JWST}
  if(filter=='F770W_JWST'){filt_F770W_JWST=NULL; data('filt_F770W_JWST', envir = environment()); out = filt_F770W_JWST}
  if(filter=='F1000W_JWST'){filt_F1000W_JWST=NULL; data('filt_F1000W_JWST', envir = environment()); out = filt_F1000W_JWST}
  if(filter=='F1130W_JWST'){filt_F1130W_JWST=NULL; data('filt_F1130W_JWST', envir = environment()); out = filt_F1130W_JWST}
  if(filter=='F1280W_JWST'){filt_F1280W_JWST=NULL; data('filt_F1280W_JWST', envir = environment()); out = filt_F1280W_JWST}
  if(filter=='F1500W_JWST'){filt_F1500W_JWST=NULL; data('filt_F1500W_JWST', envir = environment()); out = filt_F1500W_JWST}
  if(filter=='F1800W_JWST'){filt_F1800W_JWST=NULL; data('filt_F1800W_JWST', envir = environment()); out = filt_F1800W_JWST}
  if(filter=='F2100W_JWST'){filt_F2100W_JWST=NULL; data('filt_F2100W_JWST', envir = environment()); out = filt_F2100W_JWST}
  if(filter=='F2550W_JWST'){filt_F2550W_JWST=NULL; data('filt_F2550W_JWST', envir = environment()); out = filt_F2550W_JWST}

  if(filter=='R062_WFIRST'){filt_R062_WFIRST=NULL; data('filt_R062_WFIRST', envir = environment()); out = filt_R062_WFIRST}
  if(filter=='Z087_WFIRST'){filt_Z087_WFIRST=NULL; data('filt_Z087_WFIRST', envir = environment()); out = filt_Z087_WFIRST}
  if(filter=='Y106_WFIRST'){filt_Y106_WFIRST=NULL; data('filt_Y106_WFIRST', envir = environment()); out = filt_Y106_WFIRST}
  if(filter=='J129_WFIRST'){filt_J129_WFIRST=NULL; data('filt_J129_WFIRST', envir = environment()); out = filt_J129_WFIRST}
  if(filter=='W146_WFIRST'){filt_W146_WFIRST=NULL; data('filt_W146_WFIRST', envir = environment()); out = filt_W146_WFIRST}
  if(filter=='H158_WFIRST'){filt_H158_WFIRST=NULL; data('filt_H158_WFIRST', envir = environment()); out = filt_H158_WFIRST}
  if(filter=='F184_WFIRST'){filt_F184_WFIRST=NULL; data('filt_F184_WFIRST', envir = environment()); out = filt_F184_WFIRST}

  if(filter=='VIS_Euclid'){filt_VIS_Euclid=NULL; data('filt_VIS_Euclid', envir = environment()); out = filt_VIS_Euclid}
  if(filter=='Y_Euclid'){filt_Y_Euclid=NULL; data('filt_Y_Euclid', envir = environment()); out = filt_Y_Euclid}
  if(filter=='Blue_Euclid'){filt_Blue_Euclid=NULL; data('filt_Blue_Euclid', envir = environment()); out = filt_Blue_Euclid}
  if(filter=='J_Euclid'){filt_J_Euclid=NULL; data('filt_J_Euclid', envir = environment()); out = filt_J_Euclid}
  if(filter=='Red_Euclid'){filt_Red_Euclid=NULL; data('filt_Red_Euclid', envir = environment()); out = filt_Red_Euclid}
  if(filter=='H_Euclid'){filt_H_Euclid=NULL; data('filt_H_Euclid', envir = environment()); out = filt_H_Euclid}

  if(filter=='W1_WISE' | filter=='W1'){filt_W1_WISE=NULL; data('filt_W1_WISE', envir = environment()); out = filt_W1_WISE}
  if(filter=='W2_WISE' | filter=='W2'){filt_W2_WISE=NULL; data('filt_W2_WISE', envir = environment()); out = filt_W2_WISE}
  if(filter=='W3_WISE' | filter=='W3'){filt_W3_WISE=NULL; data('filt_W3_WISE', envir = environment()); out = filt_W3_WISE}
  if(filter=='W4_WISE' | filter=='W4'){filt_W4_WISE=NULL; data('filt_W4_WISE', envir = environment()); out = filt_W4_WISE}

  if(filter=='I1_Spitzer' | filter=='I1'){filt_I1_Spitzer=NULL; data('filt_I1_Spitzer', envir = environment()); out = filt_I1_Spitzer}
  if(filter=='I2_Spitzer' | filter=='I2'){filt_I2_Spitzer=NULL; data('filt_I2_Spitzer', envir = environment()); out = filt_I2_Spitzer}
  if(filter=='I3_Spitzer' | filter=='I3'){filt_I3_Spitzer=NULL; data('filt_I3_Spitzer', envir = environment()); out = filt_I3_Spitzer}
  if(filter=='I4_Spitzer' | filter=='I4'){filt_I4_Spitzer=NULL; data('filt_I4_Spitzer', envir = environment()); out = filt_I4_Spitzer}
  if(filter=='M24_Spitzer' | filter=='M24'){filt_M24_Spitzer=NULL; data('filt_M24_Spitzer', envir = environment()); out = filt_M24_Spitzer}
  if(filter=='M70_Spitzer' | filter=='M70'){filt_M70_Spitzer=NULL; data('filt_M70_Spitzer', envir = environment()); out = filt_M70_Spitzer}
  if(filter=='M160_Spitzer' | filter=='M160'){filt_M160_Spitzer=NULL; data('filt_M160_Spitzer', envir = environment()); out = filt_M160_Spitzer}

  if(filter=='P70_Herschel' | filter=='P70'){filt_P70_Herschel=NULL; data('filt_P70_Herschel', envir = environment()); out = filt_P70_Herschel}
  if(filter=='P100_Herschel' | filter=='P100'){filt_P100_Herschel=NULL; data('filt_P100_Herschel', envir = environment()); out = filt_P100_Herschel}
  if(filter=='P160_Herschel' | filter=='P160'){filt_P160_Herschel=NULL; data('filt_P160_Herschel', envir = environment()); out = filt_P160_Herschel}
  if(filter=='S250_Herschel' | filter=='S250'){filt_S250_Herschel=NULL; data('filt_S250_Herschel', envir = environment()); out = filt_S250_Herschel}
  if(filter=='S350_Herschel' | filter=='S350'){filt_S350_Herschel=NULL; data('filt_S350_Herschel', envir = environment()); out = filt_S350_Herschel}
  if(filter=='S500_Herschel' | filter=='S500'){filt_S500_Herschel=NULL; data('filt_S500_Herschel', envir = environment()); out = filt_S500_Herschel}

  if(filter=='S450_JCMT' | filter=='S450'){filt_S450_JCMT=NULL; data('filt_S450_JCMT', envir = environment()); out = filt_S450_JCMT}
  if(filter=='S850_JCMT' | filter=='S850'){filt_S850_JCMT=NULL; data('filt_S850_JCMT', envir = environment()); out = filt_S850_JCMT}

  if(filter=='1mm_Aztec' | filter=='1mm'){filt_1mm_Aztec=NULL; data('filt_1mm_Aztec', envir = environment()); out = filt_1mm_Aztec}
  if(filter=='2mm_Gizmo' | filter=='2mm'){filt_2mm_Gizmo=NULL; data('filt_2mm_Gizmo', envir = environment()); out = filt_2mm_Gizmo}

  if(filter=='Band9_ALMA' | filter=='Band9'){out = .tophat(4000000.0,5000000.0)}
  if(filter=='Band8_ALMA' | filter=='Band8'){out = .tophat(6000000.0,8000000.0)}
  if(filter=='Band7_ALMA' | filter=='Band7'){out = .tophat(8000000.0,11000000.0)}
  if(filter=='Band6_ALMA' | filter=='Band6'){out = .tophat(11000000.0,14000000.0)}
  if(filter=='Band5_ALMA' | filter=='Band5'){out = .tophat(14000000.0,18000000.0)}
  if(filter=='Band4_ALMA' | filter=='Band4'){out = .tophat(18000000.0, 24000000.0)}

  if(length(grep('tophat',filter)) == 1){
    wavelims = as.numeric(strsplit(filter,'_')[[1]][2:3])
    out = .tophat(wavelims[1], wavelims[2])
  }

  if(is.null(out)){
    EAZY_filters=NULL
    data('EAZY_filters', envir = environment())
    Eazy_try = grep(filter, EAZY_filters$info)
    if(length(Eazy_try)==1){
      out = EAZY_filters$filters[[Eazy_try]]
    }else if(length(Eazy_try)>1){
      message(paste(c('Eazy filter name is ambiguous, currently grep-ing:', EAZY_filters$info[Eazy_try]), collapse='\n'))
    }
  }
  if(is.null(out)){
    message(paste('Filter name',filter,'not recognised!'))
  }
  return(out)
}
