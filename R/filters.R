.tophat = function(wavelo, wavehi){
  return(cbind(wave=c(wavelo,wavehi), response=c(1,1)))
}

getfilt=function(filter){
  out = NULL

  valid_filters = c(
    # GALEX
    "FUV_GALEX", "NUV_GALEX",

    # SDSS
    "u_SDSS", "g_SDSS", "r_SDSS", "i_SDSS", "z_SDSS",

    # VST
    "u_VST", "g_VST", "r_VST", "i_VST", "z_VST",

    # LSST / Rubin
    "u_LSST", "g_LSST", "r_LSST", "i_LSST", "z_LSST", "y_LSST",

    # HSC
    "g_HSC", "r_HSC", "i_HSC", "z_HSC", "Y_HSC",

    # VISTA
    "Z_VISTA", "Y_VISTA", "J_VISTA", "H_VISTA", "K_VISTA", "Ks_VISTA",

    # UKIRT
    "Z_UKIRT", "Y_UKIRT", "J_UKIRT", "H_UKIRT", "K_UKIRT",

    # JWST NIRCam/MIRI
    "F070W_JWST", "F090W_JWST", "F115W_JWST", "F140M_JWST", "F150W_JWST",
    "F162M_JWST", "F164N_JWST", "F150W2_JWST", "F182M_JWST", "F187N_JWST",
    "F200W_JWST", "F210M_JWST", "F212N_JWST", "F250M_JWST", "F277W_JWST",
    "F300M_JWST", "F323N_JWST", "F322W2_JWST", "F335M_JWST", "F356W_JWST",
    "F360M_JWST", "F405N_JWST", "F410M_JWST", "F430M_JWST", "F444W_JWST",
    "F460M_JWST", "F466N_JWST", "F470N_JWST", "F480M_JWST", "F560W_JWST",
    "F770W_JWST", "F1000W_JWST", "F1130W_JWST", "F1280W_JWST", "F1500W_JWST",
    "F1800W_JWST", "F2100W_JWST", "F2550W_JWST",

    # WFIRST / Roman
    "R062_WFIRST", "Z087_WFIRST", "Y106_WFIRST", "J129_WFIRST", "W146_WFIRST",
    "H158_WFIRST", "F184_WFIRST",

    # Euclid
    "VIS_Euclid", "Y_Euclid", "Blue_Euclid", "J_Euclid", "Red_Euclid", "H_Euclid",

    # WISE
    "W1_WISE", "W2_WISE", "W3_WISE", "W4_WISE",

    # Spitzer
    "I1_Spitzer", "I2_Spitzer", "I3_Spitzer", "I4_Spitzer", "M24_Spitzer",
    "M70_Spitzer", "M160_Spitzer",

    # Herschel
    "P70_Herschel", "P100_Herschel", "P160_Herschel", "S250_Herschel",
    "S350_Herschel", "S500_Herschel",

    # JCMT
    "S450_JCMT", "S850_JCMT",

    # mm
    "1mm_Aztec", "2mm_Gizmo",

    # ALMA
    "Band9_ALMA", "Band8_ALMA", "Band7_ALMA", "Band6_ALMA", "Band5_ALMA",
    "Band4_ALMA", "Band3_ALMA", "Band2_ALMA", "Band1_ALMA",

    # VLA
    "BandQ_VLA", "BandQ", "BandKa_VLA", "BandKa", "BandK_VLA", "BandK", "BandKu_VLA",
    "BandKu", "BandX_VLA", "BandX", "BandC_VLA", "BandC", "BandS_VLA", "BandS",
    "BandL_VLA", "BandL", "BandP_VLA", "BandP",

    #GMRT
    "Band5_GMRT", "Band4_GMRT", "Band3_GMRT", "Band2_GMRT",

    #ROMAN
    "F062_WFI_ROMAN", "F087_WFI_ROMAN", "F106_WFI_ROMAN", "F129_WFI_ROMAN", "F146_WFI_ROMAN", "F158_WFI_ROMAN", "F184_WFI_ROMAN", "F213_WFI_ROMAN", "W146_WFI_ROMAN",

    #PLANCK
    "PLANCK_HFI_100GHz", "PLANCK_HFI_143GHz", "PLANCK_HFI_217GHz", "PLANCK_HFI_353GHz", "PLANCK_HFI_545GHz", "PLANCK_HFI_857GHz",

    #COBE/DIRBE
    "COBE_DIRBE_60m", "COBE_DIRBE_100m", "COBE_DIRBE_140m", "COBE_DIRBE_240m"

  )

  if(filter %in% c('FUV', 'NUV')){
    filter = paste(filter,'GALEX',sep='_')
  }

  if(filter %in% c('u', 'g', 'r', 'i', 'z')){
    filter = paste(filter,'SDSS',sep='_')
  }

  if(filter %in% c('Z', 'Y', 'J', 'H', 'K', 'Ks')){
    filter = paste(filter,'VISTA',sep='_')
  }

  if(filter %in% c('W1', 'W2', 'W3', 'W4')){
    filter = paste(filter,'WISE',sep='_')
  }

  if(filter %in% c('I1', 'I2', 'I3', 'I4', 'M24', 'M70', 'M160')){
    filter = paste(filter,'Spitzer',sep='_')
  }

  if(filter %in% c('P70', 'P100', 'P160', 'S250', 'S350', 'S500')){
    filter = paste(filter,'Herschel',sep='_')
  }

  if(filter %in% c('S450', 'S850')){
    filter = paste(filter,'JCMT',sep='_')
  }

  #Rubin is the renamed LSST, so Rubin names are aliases of the LSST curves
  if(filter %in% c('u_Rubin', 'g_Rubin', 'r_Rubin', 'i_Rubin', 'z_Rubin', 'y_Rubin', 'Y_Rubin')){
    filter = sub('_Rubin$', '_LSST', filter)
  }

  #Each survey has only one y band, so y_/Y_ are the same curve within a survey
  if(filter %in% c('Y_LSST', 'y_LSST')){
    filter = 'y_LSST'
  }
  if(filter %in% c('y_HSC', 'Y_HSC')){
    filter = 'Y_HSC'
  }

  if(filter %in% c('Band1', 'Band2', 'Band3', 'Band4', 'Band5', 'Band6', 'Band7', 'Band8', 'Band9')){
    filter = paste(filter,'ALMA',sep='_')
  }

  if(filter %in% c('BandQ', 'BandKa', 'BandK', 'BandKu', 'BandX', 'BandC', 'BandS', 'BandL', 'BandP')){
    filter = paste(filter,'VLA',sep='_')
  }

  if (filter %in% valid_filters) {
    varname <- paste0("filt_", filter)
    data(list = varname, envir = environment())
    out <- get(varname, envir = environment())
  }

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
    }else{
      message(paste('Filter name',filter,'not recognised!'))
    }
  }

  return(out)
}
