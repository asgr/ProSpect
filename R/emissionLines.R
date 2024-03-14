emissionLines=function(Ha_lum=NULL, Hb_lum=NULL, Hlines_lum=NULL, All_lum=NULL,
                       SFR=NULL, Z=0.02, q=NULL, veldisp=50, LSF=NULL, z_LSF=0,
                       lumscale=21612724, log=TRUE, range=5, res=0.5, LKL10=NULL){
  if(is.null(LKL10)){
    if(!is.null(Ha_lum) | !is.null(SFR)){
      LKL10_NormHalpha=NULL
      data('LKL10_NormHalpha', envir = environment())
      LKL10=LKL10_NormHalpha
    }else if(!is.null(Hb_lum)){
      LKL10_NormHbeta=NULL
      data('LKL10_NormHbeta', envir = environment())
      LKL10=LKL10_NormHbeta
    }else if(!is.null(Hlines_lum)){
      LKL10_NormHlines=NULL
      data('LKL10_NormHlines', envir = environment())
      LKL10=LKL10_NormHlines
    }else if(!is.null(All_lum)){
      LKL10_NormAll=NULL
      data('LKL10_NormAll', envir = environment())
      LKL10=LKL10_NormAll
    }else{
      stop("Need one of the scalings to be non-NULL!")
    }
  }
  
  Zweights = interp_param(Z, LKL10$Z, log=log)
  if(is.null(q)){
    q = Z2q(Z)
  }
  qweights = interp_param(q, LKL10$q, log=log)
  
  wavegrid = .makewavegrid(LKL10$Mappings[[1]]$wave, veldisp=veldisp, range=range, res=res)
  wavecen = LKL10$Mappings[[1]]$wave
  
  if(!is.null(LSF)){
    if(is.function(LSF)){
      vel_LSF = LSF(wavecen*(1 + z_LSF)) #to get LSF dispersion in km/s into z frame
    }else if(is.matrix(LSF) | is.data.frame(LSF)){
      vel_LSF = approx(x=log10(LSF[,1]), y=LSF[,2], xout=log10(wavecen*(1 + z_LSF)), rule=2)$y
    }else{
      vel_LSF = rep(LSF, length(wavecen))
    }
  }else{
    vel_LSF = NULL
  }
  
  lum = rep(0, length(wavegrid))
  
  if(Zweights[,'wt_lo']>0 & qweights[,'wt_lo']>0){
    linelums=LKL10$Mappings[[Zweights[,'ID_lo']]][,qweights[,'ID_lo']+3]*Zweights[,'wt_lo']*qweights[,'wt_lo']
    for(i in 1:dim(LKL10$Mappings[[1]])[1]){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, LSF=vel_LSF[i], linelum=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[,'wt_lo']>0 & qweights[,'wt_hi']>0){
    linelums=LKL10$Mappings[[Zweights[,'ID_lo']]][,qweights[,'ID_hi']+3]*Zweights[,'wt_lo']*qweights[,'wt_hi']
    for(i in 1:dim(LKL10$Mappings[[1]])[1]){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, LSF=vel_LSF[i], linelum=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[,'wt_hi']>0 & qweights[,'wt_lo']>0){
    linelums=LKL10$Mappings[[Zweights[,'ID_hi']]][,qweights[,'ID_lo']+3]*Zweights[,'wt_hi']*qweights[,'wt_lo']
    for(i in 1:dim(LKL10$Mappings[[1]])[1]){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, LSF=vel_LSF[i], linelum=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[,'wt_hi']>0 & qweights[,'wt_hi']>0){
    linelums=LKL10$Mappings[[Zweights[,'ID_hi']]][,qweights[,'ID_hi']+3]*Zweights[,'wt_hi']*qweights[,'wt_hi']
    for(i in 1:dim(LKL10$Mappings[[1]])[1]){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, LSF=vel_LSF[i], linelum=linelums[i], wavegrid=wavegrid)}
  }
  
  if(!is.null(SFR)){
    Ha_lum=SFR2Lum(SFR, lumscale=lumscale)
    lum=lum*Ha_lum
  }else if(!is.null(Ha_lum)){
    lum=lum*Ha_lum
  }else if(!is.null(Hb_lum)){
    lum=lum*Hb_lum
  }else if(!is.null(Hlines_lum)){
    lum=lum*Hlines_lum
  }else if(!is.null(All_lum)){
    lum=lum*All_lum
  }else{
    stop("Need one of the scalings to be non-NULL!")
  }
  return(data.frame(wave=wavegrid, lum=lum))
}

Z2q=function(Zgas = 0.02, q0 = 2.8e7, g0 = -1.3, Z0 = 0.012){q0*(Zgas/Z0)^g0}

SFR2Lum=function(SFR=1, lumscale=21612724){
  return(SFR*lumscale)
}

.emissionLine=function(wavecen=6562.80 , veldisp=50, LSF=NULL, linelum=1, wavegrid=NULL, range=5, res=0.5){
  #ASGR: 09-02-2024 made some changes (using sel). Need to check this works!
  if(!is.null(LSF)){
    veldisp = sqrt(veldisp^2 + LSF^2) #add in quad
  }
  
  veldisp = veldisp/(.c_to_mps/1000)
  wave_sigma = wavecen*veldisp
  if(is.null(wavegrid)){
    wavegrid = wave_sigma*rep(seq(-range,range,by=res))+wavecen
  }
  ouput = numeric(length(wavegrid))
  #this is new
  sel = which(wavegrid > wavecen - range*wave_sigma & wavegrid < wavecen + range*wave_sigma)
  ouput[sel] = dnorm(wavegrid[sel], mean=wavecen, sd=wave_sigma)*linelum
  return(ouput)
}

.makewavegrid=function(wavecens, veldisp=50, range=5, res=0.5){
  #create high density regions around lines
  veldisp = veldisp/(.c_to_mps/1000)
  grid = seq(-range,range,by=res)
  temp = outer(wavecens*veldisp, grid) + wavecens
  return(sort(unique(temp)))
}

#this is too slow- need to evaluate a fixed wave grid on the fly.

.mergespec=function(base, add){
  wave1=log10(base[,1])
  wave2=log10(add[,1])
  flux1=log10(base[,2])
  flux2=log10(add[,2])
  waveout=unique(sort(c(wave1,wave2)))

  flux1=10^approxfun(wave1, flux1, rule=2, yleft=-Inf, yright=-Inf)(waveout)
  flux2=10^approxfun(wave2, flux2, rule=2, yleft=-Inf, yright=-Inf)(waveout)

  flux1[is.na(flux1)]=0
  flux2[is.na(flux2)]=0
  return(invisible(data.frame(wave=10^waveout, flux=flux1+flux2)))
}
