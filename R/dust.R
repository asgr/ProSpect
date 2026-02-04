CF=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return(exp(-tau*(wave/pivot)^pow))
}

CF_birth=function(wave, tau=1.0, pow=-0.7, pivot=5500){
  return(exp(-tau*(wave/pivot)^pow))
}

CF_screen=function(wave, tau=0.3, pow=-0.7, pivot=5500, Eb=0, L0=2175.8, LFWHM=470){
  if(Eb>0){
    return(exp(-tau*((wave/pivot)^pow + .drude(wave, Eb=Eb, L0=L0, LFWHM=LFWHM))))
  }else{
    return(exp(-tau*(wave/pivot)^pow))
  }
}

CF_birth_atten=function(wave, flux, tau=1.0, pow=-0.7, pivot=5500){
  flux_atten=CF_birth(wave, tau=tau, pow=pow, pivot=pivot)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return(list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten))
}

CF_screen_atten=function(wave, flux, tau=0.3, pow=-0.7, pivot=5500, Eb=0, L0=2175.8, LFWHM=470){
  flux_atten=CF_screen(wave, tau=tau, pow=pow, pivot=pivot, Eb=Eb, L0=L0, LFWHM=LFWHM)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return(list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten))
}

CF_atten=function(wave, flux, tau=0.3, pow=-0.7, pivot=5500){
  flux_atten=CF_screen(wave, tau=tau, pow=pow, pivot=pivot)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return(list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten))
}

.k_lambda=function(wave, beta=1.5){
  return(wave^(-beta)/(850e4^(-beta)))
}

atten_emit=function(wave, flux, tau=0.3, pow=-0.7, alpha_SF=1.5, Dale=NULL, Dale_M2L_func=NULL, waveout=NULL, Eb=0, L0=2175.8, LFWHM=470){
  atten=CF_screen_atten(wave=wave, flux=flux, tau=tau, pow=pow, Eb=Eb, L0=L0, LFWHM=LFWHM)
  emit=Dale_interp(alpha_SF=alpha_SF, AGNfrac = 0, Dale=Dale)
  emit$Aspec=emit$Aspec*atten$total_atten
  final=addspec(wave1=wave, flux1=atten$flux, wave2=emit$Wave, flux2=emit$Aspec, extrap=0, waveout=waveout)
  if(!is.null(Dale_M2L_func)){
    dustmass=atten$total_atten/Dale_M2L_func(alpha_SF)
  }else{
    dustmass=NULL
  }
  return(list(final=final, unatten=data.frame(wave=wave, flux=flux), atten=data.frame(wave=wave, flux=atten$flux),
              emit=data.frame(wave=emit$Wave, flux=emit$Aspec), total_atten=atten$total_atten, dustmass=dustmass))
}

blackbody=function(wave, Temp = 50, k850=0.077){
  A = 4*pi*.msol_to_kg*k850/.lsol_to_w/1e10
  return(A*cosplanckLawRadWave(wave/1e10, Temp=Temp))
}

blackbody_norm=function(wave, Temp = 50, z=0, norm=1){
  output=rep(0,length(wave))
  if(length(Temp)>1){
    if(length(norm)==1){norm=rep(norm,length(Temp))}
  }
  for(i in 1:length(Temp)){
    lims = cosplanckPeakWave(Temp=Temp[i])*c(1e-3,1e3)*1e10
    scale = integrate(blackbody, lims[1], lims[2], Temp=Temp[i])$value
    #need the (1 + z) so norm works when band stretching
    output = output+blackbody(wave=wave/(1 + z), Temp=Temp[i])*norm[i]/scale/(1 + z)
  }
  return(output)
}

greybody=function(wave, Temp=50, beta=1.5, k850=0.077){
  return(.k_lambda(wave=wave, beta=beta)*blackbody(wave=wave, Temp=Temp, k850=k850))
}

greybody_norm=function(wave, Temp = 50, beta=1.5, z=0, norm=1){
  output=rep(0,length(wave))
  if(length(Temp)>1){
    if(length(beta)==1){beta=rep(beta,length(Temp))}
    if(length(norm)==1){norm=rep(norm,length(Temp))}
  }
  for(i in 1:length(Temp)){
    lims = cosplanckPeakWave(Temp=Temp[i])*c(1e-3,1e3)*1e10
    scale = integrate(greybody, lims[1], lims[2], Temp=Temp[i], beta=beta[i])$value
    #need the (1 + z) so norm works when band stretching
    output = output + greybody(wave=wave/(1 + z), Temp=Temp[i], beta=beta[i])*norm[i]/scale/(1 + z)
  }
  return(output)
}

Dale_interp=function(alpha_SF=1.5, AGNfrac=0, type='NormTot', Dale=NULL){
  if(type=='Orig'){
    if(is.null(Dale)){
      Dale_Orig=NULL
      data('Dale_Orig', envir = environment())
      Dale=Dale_Orig
    }
  }
  if(type=='Msol'){
    if(is.null(Dale)){
      Dale_Msol=NULL
      data('Dale_Msol', envir = environment())
      Dale=Dale_Msol
    }
  }
  if(type=='NormTot'){
    if(is.null(Dale)){
      Dale_NormTot=NULL
      data('Dale_NormTot', envir = environment())
      Dale=Dale_NormTot
    }
  }
  if(type=='NormAGN'){
    Dale_NormAGN=NULL
    data('Dale_NormAGN', envir = environment())
    Dale=Dale_NormAGN
  }
  if(type=='NormSFR'){
    if(is.null(Dale)){
      Dale_NormSFR=NULL
      data('Dale_NormSFR', envir = environment())
      Dale=Dale_NormSFR
    }
  }

  if(AGNfrac>0 & AGNfrac<1){
    AGNinterp=interp_quick(x=AGNfrac, Dale$AGNfrac)
  }
  if(AGNfrac<1){
    SFinterp=interp_quick(x=alpha_SF, Dale$alpha_SF)
  }

  if(AGNfrac==0){
    output=rep(0,1496)
    output=output+Dale$Aspec[[1]][SFinterp['ID_lo'],]*SFinterp['wt_lo']
    output=output+Dale$Aspec[[1]][SFinterp['ID_hi'],]*SFinterp['wt_hi']
    return(invisible(data.frame(Wave=Dale$Wave, Aspec=output)))
  }
  if(AGNfrac>0 & AGNfrac<1){
    output=rep(0,1496)
    output=output+Dale$Aspec[[AGNinterp['ID_lo']]][SFinterp['ID_lo'],]*AGNinterp['wt_lo']*SFinterp['wt_lo']
    output=output+Dale$Aspec[[AGNinterp['ID_lo']]][SFinterp['ID_hi'],]*AGNinterp['wt_lo']*SFinterp['wt_hi']
    output=output+Dale$Aspec[[AGNinterp['ID_hi']]][SFinterp['ID_lo'],]*AGNinterp['wt_hi']*SFinterp['wt_lo']
    output=output+Dale$Aspec[[AGNinterp['ID_hi']]][SFinterp['ID_hi'],]*AGNinterp['wt_hi']*SFinterp['wt_hi']
    return(invisible(data.frame(Wave=Dale$Wave, Aspec=output)))
  }
  if(AGNfrac==1){
    return(invisible(data.frame(Wave=Dale$Wave, Aspec=Dale$Aspec[[21]][1,])))
  }
}

Dale_scale=function(alpha_SF=1.5, AGNfrac=0.5, Dale_in){
  if(missing(Dale_in)){
    Dale_NormTot=NULL
    data('Dale_NormTot', envir = environment())
    Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, Dale=Dale_NormTot)
  }
  tempapproxSF=approxfun(Dale_in$Wave/1e4, Dale_in$Aspec)
  tempSFint=integrate(tempapproxSF, lower=5, upper=20)$value

  tempAGNint=3.39296e-05 #This is always the same, by definition

  AGNscale=tempSFint/(tempAGNint+tempSFint)

  NormScale=(AGNfrac*AGNscale+(1-AGNfrac)*(1-AGNscale))
  AGNfrac=(AGNfrac*AGNscale)/NormScale
  return(c(Dustfrac_bol=1-AGNfrac, AGNfrac_bol=AGNfrac))
}

dustmass=function(wave_star, lum_star_nodust, lum_star_dust, wave_dust, lum_dust){
  DustLum=sum(c(0,diff(wave_star))*(lum_star_nodust-lum_star_dust))
  #total_atten=sum(c(0,diff(wave_star))*(flux_star_nodust-flux_star_dust))
  #DustLum=total_atten/Lum2FluxFactor(z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL)
  LtoM=sum(c(0,diff(wave_dust))*lum_dust, na.rm=TRUE)
  DustMass=DustLum/LtoM
  return=c(DustMass=DustMass, DustLum=DustLum, M2L=1/LtoM)
}

.drude=function(wave, Eb=3.3, L0=2175.8, LFWHM=470){
  return(Eb*(LFWHM*wave)^2 / ((wave^2 - L0^2)^2 + (LFWHM*wave)^2))
}

Dale_M2L_variableDTH_func = function(alpha_SF, qPAH_VSG = 0.14, pivwave = 10^5.4, step_speed = 9){
  Msol = 1.989e30
  mH = 1.674e-27
  Lsol = 3.828e26
  DTH = 0.0073
  
  qBIG = 1 - qPAH_VSG
  
  ## Similar step function to what was obtained from comparing MAGPHYS mass contributing SED to total 
  smoothstep =  0.5 * (tanh(step_speed*(log10(ProSpectData::Dale_Orig$Wave) - log10(pivwave))) + 1.0)
  weight = smoothstep * (qBIG - qPAH_VSG) + qPAH_VSG
  
  ## Loop through the Dale templates and recalculate the M2L, similar to Dale_M2L
  new_M2L = sapply(
    1:64, 
    function(i){
      temp = ( (ProSpectData::Dale_Orig$Aspec[[1]][i, ] / Lsol) / (weight * DTH * mH/Msol) ) / ProSpectData::Dale_Orig$Wave
      sum( c(0, diff(ProSpectData::Dale_Orig$Wave)) * temp )
    }
  )
  
  yy = approx(
    x = ProSpectData::Dale_Orig$alpha_SF, 
    y = new_M2L, 
    xout = alpha_SF, 
    rule = 2
  )$y
  
  return( yy )
}
