SFHp4like=function(parm=c(8,9,10,10,1,0.3,1.5,0), Data, massfit=c('burstmass', 'youngmass', 'oldmass', 'ancientmass'), taufit=c('tau_birth', 'tau_screen'), dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1' , 'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose=TRUE, like=TRUE, speclib=NULL, Dale=NULL, filtout=NULL, sparse=1){
  if(zfit){
    names(parm)=c('z', massfit, taufit, dustfit)
  }else{
    names(parm)=c(massfit, taufit, dustfit)
  }
  if(verbose){print(parm)}
  if('z' %in% names(parm)){
    redshift=10^parm[names(parm)=='z']
  }else{
    redshift=Data$z
  }
  if('burstmass' %in% names(parm)){
    burstmass=10^parm[names(parm)=='burstmass']
  }else{
    burstmass = 1e+08
  }
  if('youngmass' %in% names(parm)){
    youngmass=10^parm[names(parm)=='youngmass']
  }else{
    youngmass = 1e+09
  }
  if('oldmass' %in% names(parm)){
    oldmass=10^parm[names(parm)=='oldmass']
  }else{
    oldmass = 1e+10
  }
  if('ancientmass' %in% names(parm)){
    ancientmass=10^parm[names(parm)=='ancientmass']
  }else{
    ancientmass = 1e+10
  }
  if('tau_birth' %in% names(parm)){
    tau_birth = parm[names(parm)=='tau_birth']
  }else{
    tau_birth = 1
  }
  if('tau_screen' %in% names(parm)){
    tau_screen = parm[names(parm)=='tau_screen']
  }else{
    tau_screen = 0.2
  }
  if('alpha_SF' %in% names(parm)){
    alpha_SF = parm[names(parm)=='alpha_SF']
  }else{
    alpha_SF = 1.5
  }
  if('AGNfrac' %in% names(parm)){
    AGNfrac=parm[names(parm)=='AGNfrac']
  }else{
    AGNfrac = 0
  }
  
  bad=1e10
  if(burstmass<0){if(verbose){print(-bad)}; return(bad)}
  if(youngmass<0){if(verbose){print(-bad)}; return(bad)}
  if(oldmass<0){if(verbose){print(-bad)}; return(bad)}
  if(ancientmass<0){if(verbose){print(-bad)}; return(bad)}
  if(tau_birth < 0 | tau_birth>5){if(verbose){print(-bad)}; return(bad)}
  if(tau_screen < 0 | tau_screen>3){if(verbose){print(-bad)}; return(bad)}
  if(alpha_SF<0.1 | alpha_SF>4){if(verbose){print(-bad)}; return(bad)}
  if(AGNfrac<0 | AGNfrac>0.95){if(verbose){print(-bad)}; return(bad)}
  
  if(is.null(Dale)){
    Dale_Msol=NULL
    data('Dale_Msol', envir = environment())
    Dale=Dale_Msol
  }
    
  filters=filters[filters %in% Data$flux$filter]
  Data$flux=Data$flux[Data$flux$filter %in% filters,]
  
  SFH_dust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse)
  SFH_nodust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse)
  Dale_interp=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Dale)
  dustout=dustmass(speclib$Wave, SFH_nodust$lum, SFH_dust$lum, Dale_Msol$Wave, Dale_interp)
  dustflux=Lum2Flux(Dale$Wave, Dale_interp*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  if(is.null(filtout)){
    photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  }else{
    photom_out={}
    fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
    for(i in 1:length(filtout)){
      photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
    }
  }
  
  likeout=sum(dnorm(x=(Data$flux$flux-photom_out)/Data$flux$fluxerr, log=TRUE), na.rm = TRUE)
  if(verbose){print(likeout)}
  if(like){
    return(-likeout)
  }else{
    SMtot=SMstarp5(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, speclib=speclib)
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

SFHp5like=function(parm=c(8,9,10,10,10,1,0.3,1.5,0), Data, massfit=c('burstmass', 'youngmass', 'midmass', 'oldmass', 'ancientmass'), taufit=c('tau_birth', 'tau_screen'), dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1' , 'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose=TRUE, like=TRUE, speclib=NULL, Dale=NULL, filtout=NULL, sparse=1){
  if(zfit){
    names(parm)=c('z', massfit, taufit, dustfit)
  }else{
    names(parm)=c(massfit, taufit, dustfit)
  }
  if(verbose){print(parm)}
  if('z' %in% names(parm)){
    redshift=10^parm[names(parm)=='z']
  }else{
    redshift=Data$z
  }
  if('burstmass' %in% names(parm)){
    burstmass=10^parm[names(parm)=='burstmass']
  }else{
    burstmass = 1e+08
  }
  if('youngmass' %in% names(parm)){
    youngmass=10^parm[names(parm)=='youngmass']
  }else{
    youngmass = 1e+09
  }
  if('midmass' %in% names(parm)){
    midmass=10^parm[names(parm)=='midmass']
  }else{
    midmass = 1e+10
  }
  if('oldmass' %in% names(parm)){
    oldmass=10^parm[names(parm)=='oldmass']
  }else{
    oldmass = 1e+10
  }
  if('ancientmass' %in% names(parm)){
    ancientmass=10^parm[names(parm)=='ancientmass']
  }else{
    ancientmass = 1e+10
  }
  if('tau_birth' %in% names(parm)){
    tau_birth = parm[names(parm)=='tau_birth']
  }else{
    tau_birth = 1
  }
  if('tau_screen' %in% names(parm)){
    tau_screen = parm[names(parm)=='tau_screen']
  }else{
    tau_screen = 0.2
  }
  if('alpha_SF' %in% names(parm)){
    alpha_SF = parm[names(parm)=='alpha_SF']
  }else{
    alpha_SF = 1.5
  }
  if('AGNfrac' %in% names(parm)){
    AGNfrac=parm[names(parm)=='AGNfrac']
  }else{
    AGNfrac = 0
  }
  
  bad=1e10
  if(burstmass<0){if(verbose){print(-bad)}; return(bad)}
  if(youngmass<0){if(verbose){print(-bad)}; return(bad)}
  if(midmass<0){if(verbose){print(-bad)}; return(bad)}
  if(oldmass<0){if(verbose){print(-bad)}; return(bad)}
  if(ancientmass<0){if(verbose){print(-bad)}; return(bad)}
  if(tau_birth < 0 | tau_birth>5){if(verbose){print(-bad)}; return(bad)}
  if(tau_screen < 0 | tau_screen>3){if(verbose){print(-bad)}; return(bad)}
  if(alpha_SF<0.1 | alpha_SF>4){if(verbose){print(-bad)}; return(bad)}
  if(AGNfrac<0 | AGNfrac>0.95){if(verbose){print(-bad)}; return(bad)}
  
  if(is.null(Dale)){
    Dale_Msol=NULL
    data('Dale_Msol', envir = environment())
    Dale=Dale_Msol
  }
    
  filters=filters[filters %in% Data$flux$filter]
  Data$flux=Data$flux[Data$flux$filter %in% filters,]
  
  SFH_dust=SFHp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse)
  SFH_nodust=SFHp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse)
  Dale_interp=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Dale)
  dustout=dustmass(speclib$Wave, SFH_nodust$lum, SFH_dust$lum, Dale_Msol$Wave, Dale_interp)
  dustflux=Lum2Flux(Dale$Wave, Dale_interp*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  if(is.null(filtout)){
    photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  }else{
    photom_out={}
    fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
    for(i in 1:length(filtout)){
      photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
    }
  }
  
  likeout=sum(dnorm(x=(Data$flux$flux-photom_out)/Data$flux$fluxerr, log=TRUE), na.rm = TRUE)
  if(verbose){print(likeout)}
  if(like){
    return(-likeout)
  }else{
    SMtot=SMstarp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, speclib=speclib)
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

SFHfunclike=function(parm=c(1,1,0.3,1.5,0), Data, massfunc=function(age, SFR=1){ifelse(age<1e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, massfuncfit=c('SFR'), massfuncpos=TRUE, taufit=c('tau_birth', 'tau_screen'), dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1' , 'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose=TRUE, like=TRUE, speclib=NULL, Dale=NULL, filtout=NULL, sparse=1){
  
  if(forcemass){
    if(zfit){
      names(parm)=c('z','SM', massfuncfit, taufit, dustfit)
    }else{
      names(parm)=c('SM',massfuncfit, taufit, dustfit)
    }
    forcemass=parm[names(parm)=='SM']
    parm=parm[names(parm)!='SM']
  }else{
    if(zfit){
      names(parm)=c('z', massfuncfit, taufit, dustfit)
    }else{
      names(parm)=c(massfuncfit, taufit, dustfit)
    }
  }
  
  massfuncfit=massfuncfit[massfuncfit %in% names(formals(massfunc))]
  if(length(massfuncfit)+length(taufit)+length(dustfit) != length(parm)){
    stop('massfunc inputs do not match correctly with parm!')
  }
  massfunclist=as.list(parm[1:length(massfuncfit)])
  
  if(verbose){print(parm)}
  if('z' %in% names(parm)){
    redshift=10^parm[names(parm)=='z']
  }else{
    redshift=Data$z
  }
  if('tau_birth' %in% names(parm)){
    tau_birth = parm[names(parm)=='tau_birth']
  }else{
    tau_birth = 1
  }
  if('tau_screen' %in% names(parm)){
    tau_screen = parm[names(parm)=='tau_screen']
  }else{
    tau_screen = 0.2
  }
  if('alpha_SF' %in% names(parm)){
    alpha_SF = parm[names(parm)=='alpha_SF']
  }else{
    alpha_SF = 1.5
  }
  if('AGNfrac' %in% names(parm)){
    AGNfrac=parm[names(parm)=='AGNfrac']
  }else{
    AGNfrac = 0
  }
  
  bad=1e10
  if(any(as.numeric(massfunclist)[massfuncpos]<0)){if(verbose){print(-bad)}; return(bad)}
  if(tau_birth < 0 | tau_birth>5){if(verbose){print(-bad)}; return(bad)}
  if(tau_screen < 0 | tau_screen>3){if(verbose){print(-bad)}; return(bad)}
  if(alpha_SF<0.1 | alpha_SF>4){if(verbose){print(-bad)}; return(bad)}
  if(AGNfrac<0 | AGNfrac>0.95){if(verbose){print(-bad)}; return(bad)}
  
  if(is.null(Dale)){
    Dale_Msol=NULL
    data('Dale_Msol', envir = environment())
    Dale=Dale_Msol
  }
    
  filters=filters[filters %in% Data$flux$filter]
  Data$flux=Data$flux[Data$flux$filter %in% filters,]
  
  #SFH_dust=SFHfunc(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, z=redshift, outtype=NULL, speclib=speclib)
  SFH_dust=do.call('SFHfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, tau_birth=tau_birth, tau_screen=tau_screen, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse), massfunclist))
  if(sum(SFH_dust$massvec,na.rm=TRUE)==0){if(verbose){print(-bad)}; return(bad)}
  SFH_nodust=do.call('SFHfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=speclib, sparse=sparse), massfunclist))
  Dale_interp=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Dale)
  dustout=dustmass(speclib$Wave, SFH_nodust$lum, SFH_dust$lum, Dale_Msol$Wave, Dale_interp)
  dustflux=Lum2Flux(Dale$Wave, Dale_interp*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  if(is.null(filtout)){
    photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  }else{
    photom_out={}
    fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
    for(i in 1:length(filtout)){
      photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
    }
  }
  
  likeout=sum(dnorm(x=(Data$flux$flux-photom_out)/Data$flux$fluxerr, log=TRUE), na.rm = TRUE)
  if(verbose){print(likeout)}
  if(like){
    return(-likeout)
  }else{
    SMtot=do.call('SMstarfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, z=redshift), massfunclist))
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

