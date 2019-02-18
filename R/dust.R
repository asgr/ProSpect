CF=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_birth=function(wave, tau=1.0, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_screen=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_birth_atten=function(wave, flux, tau=1.0, pow=-0.7, pivot=5500){
  flux_atten=CF_birth(wave, tau=tau, pow=pow, pivot=pivot)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return=list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten)
}

CF_screen_atten=function(wave, flux, tau=0.3, pow=-0.7, pivot=5500){
  flux_atten=CF_screen(wave, tau=tau, pow=pow, pivot=pivot)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return=list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten)
}

CF_atten=function(wave, flux, tau=0.3, pow=-0.7, pivot=5500){
  flux_atten=CF_screen(wave, tau=tau, pow=pow, pivot=pivot)*flux
  unatten=sum(flux*c(0,diff(wave)))
  atten=sum(flux_atten*c(0,diff(wave)))
  total_atten=unatten-atten
  return=list(flux=flux_atten, total_atten=total_atten, attenfrac=atten/unatten)
}

.k_lambda=function(wave, beta=1.5){
  return=wave^(-beta)/(850e4^(-beta))
}

blackbody=function(wave, Temp = 50, k850=0.077){
  A = 4*pi*.msol_to_kg*k850/.lsol_to_w/1e10
  return=A*cosplanckLawRadWave(wave/1e10, Temp=Temp)
}

blackbody_norm=function(wave, Temp = 50, z=0, norm=1){
  output=rep(0,length(wave))
  if(length(Temp)>1){
    if(length(norm)==1){norm=rep(norm,length(Temp))}
  }
  for(i in 1:length(Temp)){
    lims=cosplanckPeakWave(Temp=Temp[i])*c(1e-3,1e3)*1e10
    scale=integrate(blackbody, lims[1], lims[2], Temp=Temp[i])$value
    output=output+blackbody(wave=wave/(1+z), Temp=Temp[i])*norm[i]/scale
  }
  return=output
}

greybody=function(wave, Temp=50, beta=1.5, k850=0.077){
  return=.k_lambda(wave=wave, beta=beta)*blackbody(wave=wave, Temp=Temp, k850=k850)
}

greybody_norm=function(wave, Temp = 50, beta=1.5, z=0, norm=1){
  output=rep(0,length(wave))
  if(length(Temp)>1){
    if(length(beta)==1){beta=rep(beta,length(Temp))}
    if(length(norm)==1){norm=rep(norm,length(Temp))}
  }
  for(i in 1:length(Temp)){
    lims=cosplanckPeakWave(Temp=Temp[i])*c(1e-3,1e3)*1e10
    scale=integrate(greybody, lims[1], lims[2], Temp=Temp[i], beta=beta[i])$value
    output=output+greybody(wave=wave/(1+z), Temp=Temp[i], beta=beta[i])*norm[i]/scale
  }
  return=output
}

Dale_interp=function(alpha_SF=1.5, AGNfrac=0, type='Msol', Dale=NULL){
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
    AGNinterp=interp_param(x=AGNfrac, Dale$AGNfrac)
  }
  if(AGNfrac<1){
    SFinterp=interp_param(x=alpha_SF, Dale$alpha_SF)
  }
  
  # if(AGNfrac<0){AGNfrac=0}
  # if(AGNfrac>1){AGNfrac=1}
  # if(alpha_SF<0.0625){alpha_SF=0.0625}
  # if(alpha_SF>4.0000){alpha_SF=4.0000}
  # 
  # AGNfracloc=max(which(Dale$AGNfrac<=AGNfrac))
  # alpha_SFloc=max(which(Dale$alpha_SF<=alpha_SF))
  # 
  # if(AGNfracloc==21){AGNfracloc=20}
  # if(alpha_SFloc==64){alpha_SFloc=63}
  # 
  # AGNfraclo=Dale$AGNfrac[AGNfracloc]
  # AGNfrachi=Dale$AGNfrac[AGNfracloc+1]
  # alpha_SFlo=Dale$alpha_SF[alpha_SFloc]
  # alpha_SFhi=Dale$alpha_SF[alpha_SFloc+1]
  # 
  # AGNfraclow=(AGNfrachi-AGNfrac)/(AGNfrachi-AGNfraclo)
  # AGNfrachiw=1-AGNfraclow
  # alpha_SFlow=(alpha_SFhi-alpha_SF)/(alpha_SFhi-alpha_SFlo)
  # alpha_SFhiw=1-alpha_SFlow
  
  if(AGNfrac==0){
    output=rep(0,1496)
    output=output+Dale$Aspec[[1]][SFinterp$ID_lo,]*SFinterp$weight_lo
    output=output+Dale$Aspec[[1]][SFinterp$ID_hi,]*SFinterp$weight_hi
    return(invisible(data.frame(Wave=Dale$Wave, Aspec=output)))
  }
  if(AGNfrac>0 & AGNfrac<1){
    output=rep(0,1496)
    output=output+Dale$Aspec[[AGNinterp$ID_lo]][SFinterp$ID_lo,]*AGNinterp$weight_lo*SFinterp$weight_lo
    output=output+Dale$Aspec[[AGNinterp$ID_lo]][SFinterp$ID_hi,]*AGNinterp$weight_lo*SFinterp$weight_hi
    output=output+Dale$Aspec[[AGNinterp$ID_hi]][SFinterp$ID_lo,]*AGNinterp$weight_hi*SFinterp$weight_lo
    output=output+Dale$Aspec[[AGNinterp$ID_hi]][SFinterp$ID_hi,]*AGNinterp$weight_hi*SFinterp$weight_hi
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
