ProSpectSED=function(SFH, z=0.1, tau_birth=1, tau_screen=0.3, tau_AGN=1, pow_birth = -0.7,
                     pow_screen = -0.7, pow_AGN = -0.7, alpha_SF_birth=1, alpha_SF_screen=3,
                     alpha_SF_AGN=0, AGNlum=0, sparse=5, speclib=NULL, Dale=NULL, AGN=NULL,
                     filtout=NULL, returnall=TRUE, H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM,
                     ref, ...){
  
  tau_birth=.interval(tau_birth,0,10)
  tau_screen=.interval(tau_screen,0,10)
  tau_AGN=.interval(tau_AGN,0,10)
  pow_birth=.interval(pow_birth,0,10)
  pow_screen=.interval(pow_screen,0,10)
  pow_AGN=.interval(pow_AGN,0,10)
  alpha_SF_birth=.interval(alpha_SF_birth,0,10)
  alpha_SF_screen=.interval(alpha_SF_screen,0,10)
  alpha_SF_AGN=.interval(alpha_SF_AGN,0,10)

  Stars=SFH(z=-1, tau_birth=tau_birth, tau_screen=tau_screen, pow_birth=pow_birth, pow_screen=pow_screen, sparse=sparse, speclib=speclib, ...)
  
  Dust_Birth=Dale_interp(alpha_SF=alpha_SF_birth, Dale=Dale)
  Dust_Screen=Dale_interp(alpha_SF=alpha_SF_screen, Dale=Dale)
  
  SED_Bdust_Sdust=Dust_Screen$Aspec*Stars$lumtot_screen+Dust_Birth$Aspec*Stars$lumtot_birth
  SED_Stars_Bdust_Sdust=addspec(wave1 = Stars$wave, flux1 = Stars$lum, wave2 = Dust_Screen$Wave, flux2 = SED_Bdust_Sdust, extrap = 0)
  
  StarsAtten=data.frame(wave=Stars$wave, lum=Stars$lum)
  StarsUnAtten=data.frame(wave=Stars$wave, lum=Stars$lum_unatten)
  DustEmit=data.frame(wave=Dust_Screen$Wave, lum=SED_Bdust_Sdust)
  
  if(is.null(AGN) | AGNlum==0){
    Final=SED_Stars_Bdust_Sdust
    AGN=NULL
  }else{
    #First we attenuate by the hot taurus
    AGN=atten_emit(wave=AGN$Wave, flux=AGN$Aspec*AGNlum/(3.828e33), tau=tau_AGN, pow=pow_AGN, alpha_SF=alpha_SF_AGN, Dale=Dale)$final
    #Second we re-attenuate the above by the screen (since it still has to pass out of the galaxy)
    AGN=atten_emit(wave=AGN$wave, flux=AGN$flux, tau=tau_screen, pow=pow_screen, alpha_SF=alpha_SF_screen, Dale=Dale)$final
    Final=addspec(wave1=SED_Stars_Bdust_Sdust$wave, flux1=SED_Stars_Bdust_Sdust$flux, wave2=AGN$wave, flux2=AGN$flux, extrap=0)
    colnames(AGN)[2]='lum'
  }
  
  colnames(Final)[2]='lum'
  
  Flux=Lum2Flux(wave=Final$wave, lum=Final$lum, z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)
  
  Flux$flux=convert_wave2freq(flux_wave=Flux$flux*1e23, wave=Flux$wave)
  photom_out={}
  for(i in 1:length(filtout)){
    photom_out = c(photom_out, bandpass(flux=Flux$flux, wave=Flux$wave, filter=filtout[[i]], lum=TRUE))
  }
  
  if(returnall){
    return(invisible(list(Photom=photom_out, FinalFlux=Flux, FinalLum=Final, StarsAtten=StarsAtten, StarsUnAtten=StarsUnAtten, DustEmit=DustEmit, AGN=AGN)))
  }else{
    return(invisible(photom_out))
  }
}

ProSpectSEDlike=function(parm=c(8,9,10,10,0,-0.5,0.2), Data){
  if(Data$verbose){print(parm)}
  if(is.null(Data$logged)){
    parmlist=parm
  }else if(length(Data$logged)==1){
    if(Data$logged){
      parmlist=as.list(10^parm)
    }else{
      parmlist=parm
    }
  }else{
    parmlist=parm
    parmlist[Data$logged]=10^parmlist[Data$logged]
  }
  
  names(parmlist)=Data$parmnames
  SEDout=do.call('ProSpectSED', args=c(parmlist, list(SFH=Data$SFH), list(speclib=Data$speclib), list(Dale=Data$Dale), list(AGN=Data$AGN), list(filtout=Data$filtout), Data$fixed))
  
  cutsig=(Data$flux$flux-SEDout$Photom)/Data$flux$fluxerr
  if(Data$like=='norm'){
    LL=sum(dnorm(x=cutsig, log=TRUE), na.rm = TRUE)
  }else if(Data$like=='chisq'){
    LL=dchisq(sum(cutsig^2), df=length(Data$filtout)-length(parm), log=TRUE)
  }else if(Data$like=='st'){
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    LL=sum(dt(cutsig, df=dof, log=TRUE), na.rm = TRUE)
  }else{
    stop('Bad like option!')
  }
  if(is.null(Data$prior)){
    LP=LL
  }else{
    LP=LL+Data$prior(parm)
  }
  if(Data$verbose){print(LP)}
  if(Data$fit=='optim'){
    return(-LP)
  }
  if(Data$fit=='LD'){
    return(list(LP=LP,Dev=-2*LL,Monitor='',yhat=1,parm=parm))
  }
  if(Data$fit=='check'){
    return(invisible(list(LP=LP,Dev=-2*LL,Monitor='',yhat=1,parm=parm,SEDout=SEDout)))
  }
}