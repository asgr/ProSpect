ProSpectSED=function(SFH, z=0.1, tau_birth=1, tau_screen=0.3, tau_AGN=1, pow_birth = -0.7,
                     pow_screen = -0.7, pow_AGN = -0.7, alpha_SF_birth=1, alpha_SF_screen=3,
                     alpha_SF_AGN=0, AGNlum=0, sparse=5, speclib=NULL, Dale=NULL, AGN=NULL,
                     filtout=NULL, Dale_M2L_func=NULL, returnall=TRUE, H0=67.8,
                     OmegaM=0.308, OmegaL=1-OmegaM, waveout=seq(2,7,by=0.01), ref, unimax=13.8e9, agemax=NULL, ...){
  
  tau_birth=.interval(tau_birth,0,10,reflect=FALSE)
  tau_screen=.interval(tau_screen,0,10,reflect=FALSE)
  tau_AGN=.interval(tau_AGN,0,10,reflect=FALSE)
  pow_birth=.interval(pow_birth,-2,0,reflect=FALSE)
  pow_screen=.interval(pow_screen,-2,0,reflect=FALSE)
  pow_AGN=.interval(pow_AGN,-2,0,reflect=FALSE)
  alpha_SF_birth=.interval(alpha_SF_birth,0.0625,4,reflect=FALSE)
  alpha_SF_screen=.interval(alpha_SF_screen,0.0625,4,reflect=FALSE)
  alpha_SF_AGN=.interval(alpha_SF_AGN,0.0625,4,reflect=FALSE)

  Stars=SFH(z=z, tau_birth=tau_birth, tau_screen=tau_screen, pow_birth=pow_birth, pow_screen=pow_screen, sparse=sparse, speclib=speclib, filters=NULL, unimax=unimax, agemax=agemax, ...)
  
  Dust_Birth=Dale_interp(alpha_SF=alpha_SF_birth, Dale=Dale)
  Dust_Screen=Dale_interp(alpha_SF=alpha_SF_screen, Dale=Dale)
  
  SED_Bdust_Sdust=Dust_Birth$Aspec*Stars$lumtot_birth+Dust_Screen$Aspec*Stars$lumtot_screen
  SED_Stars_Bdust_Sdust=addspec(wave1 = Stars$wave_lum, flux1 = Stars$lum_atten, wave2 = Dust_Screen$Wave, flux2 = SED_Bdust_Sdust, extrap = 0, waveout=waveout)
  
  if(!is.null(Dale_M2L_func) & returnall){
    dustlum_birth=Stars$lumtot_birth
    dustlum_screen=Stars$lumtot_screen
    dustmass_birth=Stars$lumtot_birth/Dale_M2L_func(alpha_SF_birth)
    dustmass_screen=Stars$lumtot_screen/Dale_M2L_func(alpha_SF_screen)
  }else{
    dustlum_birth=0
    dustlum_screen=0
    dustmass_birth=0
    dustmass_screen=0
  }
  
  
  if(is.null(AGN) | AGNlum==0){
    Final=SED_Stars_Bdust_Sdust
    AGN=NULL
    dustlum_AGN=0
    dustmass_AGN=0
  }else{
    #First we attenuate by the hot taurus
    AGN=atten_emit(wave=AGN$Wave, flux=AGN$Aspec*AGNlum/(.lsun_to_erg), tau=tau_AGN, pow=pow_AGN, alpha_SF=alpha_SF_AGN, Dale=Dale, Dale_M2L_func=Dale_M2L_func, waveout=waveout)
    if(!is.null(Dale_M2L_func) & returnall){
      dustlum_AGN=AGN$total_atten
      dustmass_AGN=AGN$dustmass
    }
    #Second we re-attenuate the above by the screen (since it still has to pass out of the galaxy)
    AGN=atten_emit(wave=AGN$final$wave, flux=AGN$final$flux, tau=tau_screen, pow=pow_screen, alpha_SF=alpha_SF_screen, Dale=Dale, Dale_M2L_func=Dale_M2L_func, waveout=waveout)
    if(!is.null(Dale_M2L_func) & returnall){
      dustlum_screen=dustlum_screen+AGN$total_atten
      dustmass_screen=dustmass_screen+AGN$dustmass
    }else{
      dustlum_AGN=0
      dustmass_AGN=0
    }
    AGN=AGN$final
    Final=data.frame(wave=SED_Stars_Bdust_Sdust$wave, flux=SED_Stars_Bdust_Sdust$flux+AGN$flux)
    colnames(AGN)[2]='lum'
  }
  
  colnames(Final)[2]='lum'
  
  if(z>0 & !is.null(filtout)){
    Flux=Lum2Flux(wave=Final$wave, lum=Final$lum, z=z, H0=H0, OmegaM=OmegaM, OmegaL=OmegaL, ref=ref)
    Flux$flux=convert_wave2freq(flux_wave=Flux$flux*.cgs_to_jansky, wave=Flux$wave)
    photom_out={}
    for(i in 1:length(filtout)){
      photom_out = c(photom_out, bandpass(flux=Flux$flux, wave=Flux$wave, filter=filtout[[i]], lum=TRUE))
    }
  }else{
    Flux=NULL
    photom_out=NULL
  }
  
  if(returnall){
    StarsAtten=data.frame(wave=Stars$wave_lum, lum=Stars$lum_atten)
    StarsUnAtten=data.frame(wave=Stars$wave, lum=Stars$lum_unatten)
    DustEmit=data.frame(wave=Dust_Screen$Wave, lum=SED_Bdust_Sdust)
    return(invisible(list(Photom=photom_out, FinalFlux=Flux, FinalLum=Final, StarsAtten=StarsAtten, StarsUnAtten=StarsUnAtten, DustEmit=DustEmit, AGN=AGN, Stars=Stars, dustmass=c(birth=dustmass_birth, screen=dustmass_screen, AGN=dustmass_AGN, total=sum(c(dustmass_birth,dustmass_screen,dustmass_AGN),na.rm=TRUE)), dustlum=c(birth=dustlum_birth, screen=dustlum_screen, AGN=dustlum_AGN, total=sum(c(dustlum_birth,dustlum_screen,dustlum_AGN),na.rm=TRUE)))))
  }else{
    return(invisible(photom_out))
  }
}

ProSpectSEDlike=function(parm=c(8,9,10,10,0,-0.5,0.2), Data){
  if(is.null(Data$fit)){Data$fit='optim'}
  Data$fit=tolower(Data$fit)
  if(is.null(Data$like)){Data$like='st'}
  if(is.null(Data$mon.names)){Data$mon.names='NULL'}
  if(is.null(Data$verbose)){Data$verbose=FALSE}
  
  if(Data$fit=='check' | ((length(grep('dustmass',Data$mon.names))>0 | length(grep('dustlum',Data$mon.names))>0) & (Data$fit=='ld' | Data$fit=='la'))){
    returnall=TRUE
  }else{
    returnall=FALSE
  }
  
  names(parm)=Data$parm.names
  
  if(!is.null(Data$constraints)){
    parm=Data$constraints(parm)
  }
  
  if(!is.null(Data$intervals)){
    parm[parm<Data$intervals$lo]=Data$intervals$lo[parm<Data$intervals$lo]
    parm[parm>Data$intervals$hi]=Data$intervals$hi[parm>Data$intervals$hi]
  }
  
  if(!is.null(Data$logged)){
    if(length(Data$logged)==1){
      if(Data$logged){
        parmlist=10^parm
      }else{
        parmlist=parm
      }
    }else{
      parmlist=parm
      parmlist[Data$logged]=10^parm[Data$logged]
    }
  }else{
    parmlist=parm
  }
  
  if(Data$verbose){print(parmlist)}
  
  if(returnall){
    SEDout=do.call('ProSpectSED', args=c(parmlist, list(SFH=Data$SFH), list(speclib=Data$speclib), list(Dale=Data$Dale), list(AGN=Data$AGN), list(filtout=Data$filtout), list(returnall=TRUE), list(Dale_M2L_func=Data$Dale_M2L_func), Data$arglist))
    Monitor={}
    if(length(grep('dustmass',Data$mon.names))>0){
      Monitor=c(dustmass=SEDout$dustmass)
    }
    if(length(grep('dustlum',Data$mon.names))>0){
      Monitor=c(dustlum=SEDout$dustlum,Monitor)
    }
    if('masstot' %in% Data$mon.names){
      Monitor=c(masstot=SEDout$Stars$masstot,Monitor)
    }
    if('fluxes' %in% Data$mon.names){
      Monitor=c(SEDout$Photom,Monitor)
    }
    Photom=SEDout$Photom
  }else{
    Photom=do.call('ProSpectSED', args=c(parmlist, list(SFH=Data$SFH), list(speclib=Data$speclib), list(Dale=Data$Dale), list(AGN=Data$AGN), list(filtout=Data$filtout), list(returnall=FALSE), Data$arglist))
    Monitor=NULL
  }
  
  cutsig=(Data$flux$flux-Photom)/Data$flux$fluxerr
  if(Data$like=='norm'){
    LL=sum(dnorm(x=cutsig, log=TRUE), na.rm = TRUE)
  }else if(Data$like=='chisq'){
    LL=dchisq(sum(cutsig^2), df=length(Data$filtout)-length(parm), log=TRUE)
  }else if(Data$like=='st'){
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    LL=as.numeric(sum(dt(cutsig, df=dof, log=TRUE), na.rm = TRUE))
  }else{
    stop('Bad like option!')
  }
  
  if(is.null(Data$prior)){
    LP=LL
  }else{
    LP=LL+as.numeric(Data$prior(parm))
  }
  if(Data$verbose){
    print(LP)
  }
  
  if(Data$fit=='ld' | Data$fit=='la' | Data$fit=='check'){
    if('LP' %in% Data$mon.names){
      Monitor=c(LP=LP,Monitor)
    }
    if(length(Monitor)==0){
      Monitor=NA
    }else{
      Monitor=Monitor[match(Data$mon.names, names(Monitor))]
    }
  }
  
  # Various returns:
  
  if(Data$fit=='optim' | Data$fit=='cma'){
    return(-LP)
  }else if(Data$fit=='ld' | Data$fit=='la'){
    return(list(LP=LP,Dev=-2*LL,Monitor=Monitor,yhat=1,parm=parm))
  }else if(Data$fit=='check'){
    return(invisible(list(LP=LP,Dev=-2*LL,Monitor=Monitor,yhat=1,parm=parm,SEDout=SEDout)))
  }else{
    return('Bad fit type!')
  }
}