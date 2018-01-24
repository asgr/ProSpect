CF_birth=function(wave, tau=1.0, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_screen=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_birth_atten=function(wave, flux, tau=1.0, pow=-0.7, pivot=5500){
  flux_atten=CF_birth(wave, tau=tau, pow=pow, pivot=pivot)*flux
  total_atten=sum((flux-flux_atten)*c(0,diff(wave)))
  return=list(flux=flux_atten, total_atten=total_atten)
}

CF_screen_atten=function(wave, flux, tau=0.3, pow=-0.7, pivot=5500){
  flux_atten=CF_screen(wave, tau=tau, pow=pow, pivot=pivot)*flux
  total_atten=sum((flux-flux_atten)*c(0,diff(wave)))
  return=list(flux=flux_atten, total_atten=total_atten)
}

.k_lambda=function(wave, beta=1.5){
  return=wave^(-beta)/(850e4^(-beta))
}

blackbody=function(wave, Temp = 50, k850=0.077){
  A = 4*pi*1.989e30*k850/3.828e26/1e10
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

Dale_interp=function(AGNfrac=0, alpha_SF=1.5, type='Msol'){
  if(type=='Orig'){
    Dale_Orig=NULL
    data('Dale_Orig', envir = environment())
    temp=Dale_Orig
  }
  if(type=='Msol'){
    Dale_Msol=NULL
    data('Dale_Msol', envir = environment())
    temp=Dale_Msol
  }
  if(type=='NormTot'){
    Dale_NormTot=NULL
    data('Dale_NormTot', envir = environment())
    temp=Dale_NormTot
  }
  if(type=='NormAGN'){
    Dale_NormAGN=NULL
    data('Dale_NormAGN', envir = environment())
    temp=Dale_NormAGN
  }
  if(type=='NormSFR'){
    Dale_NormSFR=NULL
    data('Dale_NormSFR', envir = environment())
    temp=Dale_NormSFR
  }
  
  if(AGNfrac<0){AGNfrac=0}
  if(AGNfrac>1){AGNfrac=1}
  if(alpha_SF<0.0625){alpha_SF=0.0625}
  if(alpha_SF>4.0000){alpha_SF=4.0000}
  
  AGNfracloc=max(which(temp$AGNfrac<=AGNfrac))
  alpha_SFloc=max(which(temp$alpha_SF<=alpha_SF))
  
  if(AGNfracloc==21){AGNfracloc=20}
  if(alpha_SFloc==64){alpha_SFloc=63}
  
  AGNfraclo=temp$AGNfrac[AGNfracloc]
  AGNfrachi=temp$AGNfrac[AGNfracloc+1]
  alpha_SFlo=temp$alpha_SF[alpha_SFloc]
  alpha_SFhi=temp$alpha_SF[alpha_SFloc+1]
  
  AGNfraclow=(AGNfrachi-AGNfrac)/(AGNfrachi-AGNfraclo)
  AGNfrachiw=1-AGNfraclow
  alpha_SFlow=(alpha_SFhi-alpha_SF)/(alpha_SFhi-alpha_SFlo)
  alpha_SFhiw=1-alpha_SFlow
  
  output=rep(0,1496)
  
  output=output+temp$Aspec[[AGNfracloc]][alpha_SFloc,]*AGNfraclow*alpha_SFlow
  output=output+temp$Aspec[[AGNfracloc]][alpha_SFloc+1,]*AGNfraclow*alpha_SFhiw
  output=output+temp$Aspec[[AGNfracloc+1]][alpha_SFloc,]*AGNfrachiw*alpha_SFlow
  output=output+temp$Aspec[[AGNfracloc+1]][alpha_SFloc+1,]*AGNfrachiw*alpha_SFhiw
  return=output
}
