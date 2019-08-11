emissionLines=function(Ha_lum=1, SFR=NULL, Z=0.02, q=NULL, veldisp=50, lumscale=21684047, log=TRUE, method='linear', range=5, res=0.5, LKL10=NULL){
  if(is.null(LKL10)){
    LKL10=NULL
    data('LKL10', envir = environment())
  }
  
  Zweights=interp_param(Z, LKL10$Z, log=log, method = method)
  if(is.null(q)){
    q=Z2q(Z)
  }
  qweights=interp_param(q, LKL10$q, log=log, method = method)
  
  wavegrid=.makewavegrid(LKL10$Mappings[[1]]$wave, veldisp=veldisp, range=range, res=res)
  wavecen=LKL10$Mappings[[1]]$wave
  lum=rep(0, length(wavegrid))
  
  if(Zweights[1,'weight_lo']>0 & qweights[1,'weight_lo']>0){
    linelums=LKL10$Mappings[[Zweights[1,'ID_lo']]][,qweights[1,'ID_lo']+3]*Zweights[1,'weight_lo']*qweights[1,'weight_lo']
    for(i in 1:50){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, linelums=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[1,'weight_lo']>0 & qweights[1,'weight_hi']>0){
    linelums=LKL10$Mappings[[Zweights[1,'ID_lo']]][,qweights[1,'ID_hi']+3]*Zweights[1,'weight_lo']*qweights[1,'weight_hi']
    for(i in 1:50){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, linelums=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[1,'weight_hi']>0 & qweights[1,'weight_lo']>0){
    linelums=LKL10$Mappings[[Zweights[1,'ID_hi']]][,qweights[1,'ID_lo']+3]*Zweights[1,'weight_hi']*qweights[1,'weight_lo']
    for(i in 1:50){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, linelums=linelums[i], wavegrid=wavegrid)}
  }
  if(Zweights[1,'weight_hi']>0 & qweights[1,'weight_hi']>0){
    linelums=LKL10$Mappings[[Zweights[1,'ID_hi']]][,qweights[1,'ID_hi']+3]*Zweights[1,'weight_hi']*qweights[1,'weight_hi']
    for(i in 1:50){lum=lum+.emissionLine(wavecen[i], veldisp=veldisp, linelums=linelums[i], wavegrid=wavegrid)}
  }
  
  if(!is.null(SFR)){
    Ha_lum=SFR2Lum(SFR, lumscale=lumscale)
  }
  #need to fix how the emission features stack
  return(data.frame(wave=wavegrid, lum=lum*Ha_lum))
}

Z2q=function(Zgas = 0.02, q0 = 2.8e7, g0 = -1.3, Z0 = 0.012){q0*(Zgas/Z0)^g0}

SFR2Lum=function(SFR=1, lumscale=21684047){
  return(SFR*lumscale)
}

.emissionLine=function(wavecen=6562.80 , veldisp=50, linelums=1, wavegrid=NULL, range=5, res=0.5){
  veldisp=veldisp/(.c_to_mps/1000)
  wave_sigma=wavecen*veldisp
  if(is.null(wavegrid)){
    wavegrid=wave_sigma*rep(seq(-range,range,by=res))+wavecen
  }
  return(dnorm(wavegrid, mean=wavecen, sd=wave_sigma)*linelums)
}

.makewavegrid=function(wavecens, veldisp=50, range=5, res=0.2){
  veldisp=veldisp/(.c_to_mps/1000)
  grid=rep(seq(-range,range,by=res))
  temp=outer(wavecens*veldisp, grid) + wavecens
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
