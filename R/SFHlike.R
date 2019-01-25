SFHp4like=function(parm=c(8,9,10,10,0,0.5,0.2,-2), Data, massfit=c('burstmass', 'youngmass', 'oldmass', 'ancientmass'), taufit=c('tau_birth', 'tau_screen'), powfit=NULL, dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, verbose=TRUE, sparse=5){
  
  if('massfit' %in% names(Data)){massfit = Data$massfit}
  if('taufit' %in% names(Data)){taufit = Data$taufit}
  if('powfit' %in% names(Data)){powfit = Data$powfit}
  if('dustfit' %in% names(Data)){dustfit = Data$dustfit}
  if('zfit' %in% names(Data)){zfit = Data$zfit}
  if('filters' %in% names(Data)){filters = Data$filters}
  if('verbose' %in% names(Data)){verbose = Data$verbose}
  if('sparse' %in% names(Data)){sparse = Data$sparse}
  
  if(zfit){
    names(parm)=c('z', massfit, taufit, powfit, dustfit)
  }else{
    names(parm)=c(massfit, taufit, powfit, dustfit)
  }
  if(verbose){print(parm)}
  if('z' %in% names(parm)){
    redshift=10^parm[names(parm)=='z']
  }else{
    if(!is.null(Data$fixed$z)){
      redshift=Data$fixed$z
    }else{
      redshift=0
    }
  }
  if('burstmass' %in% names(parm)){
    parm[names(parm)=='burstmass']=.interval(parm[names(parm)=='burstmass'], 0, 100, reflect=FALSE)
    burstmass=10^parm[names(parm)=='burstmass']
  }else{
    if(!is.null(Data$fixed$burstmass)){
      burstmass = Data$fixed$burstmass
    }else{
      burstmass = 0
    }
  }
  if('youngmass' %in% names(parm)){
    parm[names(parm)=='youngmass']=.interval(parm[names(parm)=='youngmass'], 0, 100, reflect=FALSE)
    youngmass=10^parm[names(parm)=='youngmass']
  }else{
    if(!is.null(Data$fixed$youngmass)){
      youngmass = Data$fixed$youngmass
    }else{
      youngmass = 0
    }
  }
  if('oldmass' %in% names(parm)){
    parm[names(parm)=='oldmass']=.interval(parm[names(parm)=='oldmass'], 0, 100, reflect=FALSE)
    oldmass=10^parm[names(parm)=='oldmass']
  }else{
    if(!is.null(Data$fixed$oldmass)){
      oldmass = Data$fixed$oldmass
    }else{
      oldmass = 0
    }
  }
  if('ancientmass' %in% names(parm)){
    parm[names(parm)=='ancientmass']=.interval(parm[names(parm)=='ancientmass'], 0, 100, reflect=FALSE)
    ancientmass=10^parm[names(parm)=='ancientmass']
  }else{
    if(!is.null(Data$fixed$ancientmass)){
      ancientmass = Data$fixed$ancientmass
    }else{
      ancientmass = 0
    }
  }
  if('tau_birth' %in% names(parm)){
    parm[names(parm)=='tau_birth']=.interval(parm[names(parm)=='tau_birth'], -2, 2, reflect=FALSE)
    tau_birth = 10^parm[names(parm)=='tau_birth']
  }else{
    if(!is.null(Data$fixed$tau_birth)){
      tau_birth = Data$fixed$tau_birth
    }else{
      tau_birth = 1
    }
  }
  if('tau_screen' %in% names(parm)){
    parm[names(parm)=='tau_screen']=.interval(parm[names(parm)=='tau_screen'], -2, 2, reflect=FALSE)
    tau_screen = 10^parm[names(parm)=='tau_screen']
  }else{
    if(!is.null(Data$fixed$tau_screen)){
      tau_screen = Data$fixed$tau_screen
    }else{
      tau_screen = 0.3
    }
  }
  if('pow_birth' %in% names(parm)){
    parm[names(parm)=='pow_birth']=.interval(parm[names(parm)=='pow_birth'], -1, 0, reflect=FALSE)
    pow_birth = 10^parm[names(parm)=='pow_birth']
  }else{
    if(!is.null(Data$fixed$pow_birth)){
      pow_birth = Data$fixed$pow_birth
    }else{
      pow_birth = -0.7
    }
  }
  if('pow_screen' %in% names(parm)){
    parm[names(parm)=='pow_screen']=.interval(parm[names(parm)=='pow_screen'], -1, 0, reflect=FALSE)
    pow_screen = 10^parm[names(parm)=='pow_screen']
  }else{
    if(!is.null(Data$fixed$pow_screen)){
      pow_screen = Data$fixed$pow_screen
    }else{
      pow_screen = -0.7
    }
  }
  if('alpha_SF' %in% names(parm)){
    parm[names(parm)=='alpha_SF']=.interval(parm[names(parm)=='alpha_SF'], -1, 0.6, reflect=FALSE)
    alpha_SF = 10^parm[names(parm)=='alpha_SF']
  }else{
    if(!is.null(Data$fixed$alpha_SF)){
      alpha_SF = Data$fixed$alpha_SF
    }else{
      alpha_SF = 1.5
    }
  }
  if('AGNfrac' %in% names(parm)){
    parm[names(parm)=='AGNfrac']=.interval(parm[names(parm)=='AGNfrac'], -2, -0.023, reflect=FALSE)
    AGNfrac = 10^parm[names(parm)=='AGNfrac']
  }else{
    if(!is.null(Data$fixed$AGNfrac)){
      AGNfrac = Data$fixed$AGNfrac
    }else{
      AGNfrac = 0
    }
  }
  
  if(!is.null(Data$Z)){
    Z=Data$Z
  }else{
    Z=5
  }
  
  # if(filters[1]=='all'){
  #   cenwave=NULL
  #   data('cenwave', envir = environment())
  #   filters=cenwave$filter
  #   filters=filters[filters %in% Data$flux$filter]
  #   Data$flux=Data$flux[Data$flux$filter %in% filters,]
  # }
  
  SFH_dust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, pow_birth=pow_birth, pow_screen=pow_screen, z=redshift, outtype=NULL, speclib=Data$speclib, Z=Z, sparse=sparse)
  SFH_nodust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=Data$speclib, Z=Z, sparse=sparse)
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Data$Dale)
  dustout=dustmass(Data$speclib$Wave[seq(1,dim(Data$speclib$Zspec[[1]])[2],by=sparse)], SFH_nodust$lum, SFH_dust$lum, Dale_in$Wave, Dale_in$Aspec)
  dustflux=Lum2Flux(Dale_in$Wave, Dale_in$Aspec*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, Dale=Data$Dale)
  Dale_in[,2]=Dale_in[,2]*dustout[3]
  Fractions=Dale_scale(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale_in=Dale_in)
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  # if(is.null(filtout)){
  #   photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  # }else{
  #   photom_out={}
  #   fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  #   for(i in 1:length(filtout)){
  #     photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
  #   }
  # }
  fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  photom_out={}
  for(i in 1:length(Data$filtout)){
    photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = Data$filtout[[i]], lum = T)*1e23)
  }
  
  cutsig=(Data$flux$flux-photom_out)/Data$flux$fluxerr
  if(Data$like=='norm'){
    likeout=sum(dnorm(x=cutsig, log=TRUE), na.rm = TRUE)
  }else if(Data$like=='chisq'){
    likeout=dchisq(sum(cutsig^2), df=Data$N-length(parm), log=TRUE)
  }else if(Data$like=='st'){
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    likeout=sum(dt(cutsig, df=dof, log=TRUE), na.rm = TRUE)
    }else{
    stop('Bad like option!')
  }
  if(verbose){print(likeout)}
  if(Data$fit=='optim'){
    return(-likeout)
  }else if(Data$fit=='LD'){
    DustMass=log10(dustout[1])
    DustLum=log10(dustout[2])
    Dustfrac_bol=log10(Fractions[1])
    AGNfrac_bol=log10(Fractions[2])
    AGNLum=log10(dustout[2]*Fractions[2]/Fractions[1])
    Monitor=c(DustMass, DustLum, AGNLum, Dustfrac_bol, AGNfrac_bol)
    Monitor[!is.finite(Monitor)]=1e10
    return(list(LP=likeout,Dev=-2*likeout,Monitor=Monitor,yhat=1,parm=parm))
  }else if(Data$fit=='check'){
    SMtot=SMstarp5(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, speclib=Data$speclib)
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

SFHp5like=function(parm=c(8,9,10,10,10,0,0.5,0.2,-2), Data, massfit=c('burstmass', 'youngmass', 'midmass', 'oldmass', 'ancientmass'), taufit=c('tau_birth', 'tau_screen'), powfit=NULL, dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, verbose=TRUE, sparse=5){
  if('massfit' %in% names(Data)){massfit = Data$massfit}
  if('taufit' %in% names(Data)){taufit = Data$taufit}
  if('powfit' %in% names(Data)){powfit = Data$powfit}
  if('dustfit' %in% names(Data)){dustfit = Data$dustfit}
  if('zfit' %in% names(Data)){zfit = Data$zfit}
  if('filters' %in% names(Data)){filters = Data$filters}
  if('verbose' %in% names(Data)){verbose = Data$verbose}
  if('sparse' %in% names(Data)){sparse = Data$sparse}
  
  if(zfit){
    names(parm)=c('z', massfit, taufit, powfit, dustfit)
  }else{
    names(parm)=c(massfit, taufit, powfit, dustfit)
  }
  if(verbose){print(parm)}
  if('z' %in% names(parm)){
    redshift=10^parm[names(parm)=='z']
  }else{
    if(!is.null(Data$fixed$z)){
      redshift=Data$fixed$z
    }else{
      redshift=0
    }
  }
  if('burstmass' %in% names(parm)){
    parm[names(parm)=='burstmass']=.interval(parm[names(parm)=='burstmass'], 0, 100, reflect=FALSE)
    burstmass=10^parm[names(parm)=='burstmass']
  }else{
    if(!is.null(Data$fixed$burstmass)){
      burstmass = Data$fixed$burstmass
    }else{
      burstmass = 0
    }
  }
  if('youngmass' %in% names(parm)){
    parm[names(parm)=='youngmass']=.interval(parm[names(parm)=='youngmass'], 0, 100, reflect=FALSE)
    youngmass=10^parm[names(parm)=='youngmass']
  }else{
    if(!is.null(Data$fixed$youngmass)){
      youngmass = Data$fixed$youngmass
    }else{
      youngmass = 0
    }
  }
  if('midmass' %in% names(parm)){
    parm[names(parm)=='midmass']=.interval(parm[names(parm)=='midmass'], 0, 100, reflect=FALSE)
    midmass=10^parm[names(parm)=='midmass']
  }else{
    if(!is.null(Data$fixed$midmass)){
      midmass = Data$fixed$midmass
    }else{
      midmass = 0
    }
  }
  if('oldmass' %in% names(parm)){
    parm[names(parm)=='oldmass']=.interval(parm[names(parm)=='oldmass'], 0, 100, reflect=FALSE)
    oldmass=10^parm[names(parm)=='oldmass']
  }else{
    if(!is.null(Data$fixed$oldmass)){
      oldmass = Data$fixed$oldmass
    }else{
      oldmass = 0
    }
  }
  if('ancientmass' %in% names(parm)){
    parm[names(parm)=='ancientmass']=.interval(parm[names(parm)=='ancientmass'], 0, 100, reflect=FALSE)
    ancientmass=10^parm[names(parm)=='ancientmass']
  }else{
    if(!is.null(Data$fixed$ancientmass)){
      ancientmass = Data$fixed$ancientmass
    }else{
      ancientmass = 0
    }
  }
  if('tau_birth' %in% names(parm)){
    parm[names(parm)=='tau_birth']=.interval(parm[names(parm)=='tau_birth'], -2, 2, reflect=FALSE)
    tau_birth = 10^parm[names(parm)=='tau_birth']
  }else{
    if(!is.null(Data$fixed$tau_birth)){
      tau_birth = Data$fixed$tau_birth
    }else{
      tau_birth = 1
    }
  }
  if('tau_screen' %in% names(parm)){
    parm[names(parm)=='tau_screen']=.interval(parm[names(parm)=='tau_screen'], -2, 2, reflect=FALSE)
    tau_screen = 10^parm[names(parm)=='tau_screen']
  }else{
    if(!is.null(Data$fixed$tau_screen)){
      tau_screen = Data$fixed$tau_screen
    }else{
      tau_screen = 0.3
    }
  }
  if('pow_birth' %in% names(parm)){
    parm[names(parm)=='pow_birth']=.interval(parm[names(parm)=='pow_birth'], -1, 0, reflect=FALSE)
    pow_birth = 10^parm[names(parm)=='pow_birth']
  }else{
    if(!is.null(Data$fixed$pow_birth)){
      pow_birth = Data$fixed$pow_birth
    }else{
      pow_birth = -0.7
    }
  }
  if('pow_screen' %in% names(parm)){
    parm[names(parm)=='pow_screen']=.interval(parm[names(parm)=='pow_screen'], -1, 0, reflect=FALSE)
    pow_screen = 10^parm[names(parm)=='pow_screen']
  }else{
    if(!is.null(Data$fixed$pow_screen)){
      pow_screen = Data$fixed$pow_screen
    }else{
      pow_screen = -0.7
    }
  }
  if('alpha_SF' %in% names(parm)){
    parm[names(parm)=='alpha_SF']=.interval(parm[names(parm)=='alpha_SF'], -1, 0.6, reflect=FALSE)
    alpha_SF = 10^parm[names(parm)=='alpha_SF']
  }else{
    if(!is.null(Data$fixed$alpha_SF)){
      alpha_SF = Data$fixed$alpha_SF
    }else{
      alpha_SF = 1.5
    }
  }
  if('AGNfrac' %in% names(parm)){
    parm[names(parm)=='AGNfrac']=.interval(parm[names(parm)=='AGNfrac'], -2, -0.023, reflect=FALSE)
    AGNfrac = 10^parm[names(parm)=='AGNfrac']
  }else{
    if(!is.null(Data$fixed$AGNfrac)){
      AGNfrac = Data$fixed$AGNfrac
    }else{
      AGNfrac = 0
    }
  }
  
  if(!is.null(Data$Z)){
    Z=Data$Z
  }else{
    Z=5
  }
  
  # bad=1e10
  # if(burstmass<0){if(verbose){print(-bad)}; return(bad)}
  # if(youngmass<0){if(verbose){print(-bad)}; return(bad)}
  # if(midmass<0){if(verbose){print(-bad)}; return(bad)}
  # if(oldmass<0){if(verbose){print(-bad)}; return(bad)}
  # if(ancientmass<0){if(verbose){print(-bad)}; return(bad)}
  # if(tau_birth < 0 | tau_birth>10){if(verbose){print(-bad)}; return(bad)}
  # if(tau_screen < 0 | tau_screen>10){if(verbose){print(-bad)}; return(bad)}
  # if(alpha_SF<0.1 | alpha_SF>4){if(verbose){print(-bad)}; return(bad)}
  # if(AGNfrac<0 | AGNfrac>0.95){if(verbose){print(-bad)}; return(bad)}
  
  # if(is.null(Dale)){
  #   Dale_Msol=NULL
  #   data('Dale_Msol', envir = environment())
  #   Dale=Dale_Msol
  # }
  
  # if(filters[1]=='all'){
  #   cenwave=NULL
  #   data('cenwave', envir = environment())
  #   filters=cenwave$filter
  #   filters=filters[filters %in% Data$flux$filter]
  #   Data$flux=Data$flux[Data$flux$filter %in% filters,]
  # }
  
  SFH_dust=SFHp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, pow_birth=pow_birth, pow_screen=pow_screen, z=redshift, outtype=NULL, speclib=Data$speclib, Z=Z, sparse=sparse)
  SFH_nodust=SFHp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=Data$speclib, Z=Z, sparse=sparse)
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Data$Dale)
  dustout=dustmass(Data$speclib$Wave[seq(1,dim(Data$speclib$Zspec[[1]])[2],by=sparse)], SFH_nodust$lum, SFH_dust$lum, Dale_in$Wave, Dale_in$Aspec)
  dustflux=Lum2Flux(Dale_in$Wave, Dale_in$Aspec*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, Dale=Data$Dale)
  Dale_in[,2]=Dale_in[,2]*dustout[3]
  Fractions=Dale_scale(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale_in=Dale_in)
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  # if(is.null(filtout)){
  #   photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  # }else{
  #   photom_out={}
  #   fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  #   for(i in 1:length(filtout)){
  #     photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
  #   }
  # }
  fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  photom_out={}
  for(i in 1:length(Data$filtout)){
    photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = Data$filtout[[i]], lum = T)*1e23)
  }
  
  cutsig=(Data$flux$flux-photom_out)/Data$flux$fluxerr
  if(Data$like=='norm'){
    likeout=sum(dnorm(x=cutsig, log=TRUE), na.rm = TRUE)
  }else if(Data$like=='chisq'){
    likeout=dchisq(sum(cutsig^2), df=Data$N-length(parm), log=TRUE)
  }else if(Data$like=='st'){
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    likeout=sum(dt(cutsig, df=dof, log=TRUE), na.rm = TRUE)
    }else{
    stop('Bad like option!')
  }
  if(verbose){print(likeout)}
  if(Data$fit=='optim'){
    return(-likeout)
  }else if(Data$fit=='LD'){
    DustMass=log10(dustout[1])
    DustLum=log10(dustout[2])
    Dustfrac_bol=log10(Fractions[1])
    AGNfrac_bol=log10(Fractions[2])
    AGNLum=log10(dustout[2]*Fractions[2]/Fractions[1])
    Monitor=c(DustMass, DustLum, AGNLum, Dustfrac_bol, AGNfrac_bol)
    Monitor[!is.finite(Monitor)]=1e10
    return(list(LP=likeout,Dev=-2*likeout,Monitor=Monitor,yhat=1,parm=parm))
  }else if(Data$fit=='check'){
    SMtot=SMstarp5(burstmass=burstmass, youngmass=youngmass, midmass=midmass, oldmass=oldmass, ancientmass=ancientmass, speclib=Data$speclib)
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

SFHfunclike=function(parm=c(1,0,0.5,0.2,-2), Data, massfunc=function(age, SFR=1){ifelse(age<1e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, massfuncfit='SFR', massfuncpos=TRUE, taufit=c('tau_birth', 'tau_screen'), powfit=NULL, dustfit=c('alpha_SF', 'AGNfrac'), zfit=FALSE, verbose=TRUE, sparse=5){
  
  if('massfunc' %in% names(Data)){massfunc = Data$massfunc}
  if('forcemass' %in% names(Data)){forcemass = Data$forcemass}
  if('unimax' %in% names(Data)){unimax = Data$unimax}
  if('agescale' %in% names(Data)){agescale = Data$agescale}
  if('massfuncfit' %in% names(Data)){massfuncfit = Data$massfuncfit}
  if('massfuncpos' %in% names(Data)){massfuncpos = Data$massfuncpos}
  if('taufit' %in% names(Data)){taufit = Data$taufit}
  if('powfit' %in% names(Data)){powfit = Data$powfit}
  if('dustfit' %in% names(Data)){dustfit = Data$dustfit}
  if('zfit' %in% names(Data)){zfit = Data$zfit}
  if('filters' %in% names(Data)){filters = Data$filters}
  if('verbose' %in% names(Data)){verbose = Data$verbose}
  if('sparse' %in% names(Data)){sparse = Data$sparse}
  
  if(forcemass){
    if(zfit){
      names(parm)=c('z','SM', massfuncfit, taufit, powfit, dustfit)
    }else{
      names(parm)=c('SM',massfuncfit, taufit, powfit, dustfit)
    }
    forcemass=parm[names(parm)=='SM']
    parm=parm[names(parm)!='SM']
  }else{
    if(zfit){
      names(parm)=c('z', massfuncfit, taufit, powfit, dustfit)
    }else{
      names(parm)=c(massfuncfit, taufit, powfit, dustfit)
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
    if(!is.null(Data$fixed$z)){
      redshift=Data$fixed$z
    }else{
      redshift=0
    }
  }
  if('tau_birth' %in% names(parm)){
    parm[names(parm)=='tau_birth']=.interval(parm[names(parm)=='tau_birth'], -2, 2, reflect=FALSE)
    tau_birth = 10^parm[names(parm)=='tau_birth']
  }else{
    if(!is.null(Data$fixed$tau_birth)){
      tau_birth = Data$fixed$tau_birth
    }else{
      tau_birth = 1
    }
  }
  if('tau_screen' %in% names(parm)){
    parm[names(parm)=='tau_screen']=.interval(parm[names(parm)=='tau_screen'], -2, 2, reflect=FALSE)
    tau_screen = 10^parm[names(parm)=='tau_screen']
  }else{
    if(!is.null(Data$fixed$tau_screen)){
      tau_screen = Data$fixed$tau_screen
    }else{
      tau_screen = 0.3
    }
  }
  if('pow_birth' %in% names(parm)){
    parm[names(parm)=='pow_birth']=.interval(parm[names(parm)=='pow_birth'], -1, 0, reflect=FALSE)
    pow_birth = 10^parm[names(parm)=='pow_birth']
  }else{
    if(!is.null(Data$fixed$pow_birth)){
      pow_birth = Data$fixed$pow_birth
    }else{
      pow_birth = -0.7
    }
  }
  if('pow_screen' %in% names(parm)){
    parm[names(parm)=='pow_screen']=.interval(parm[names(parm)=='pow_screen'], -1, 0, reflect=FALSE)
    pow_screen = 10^parm[names(parm)=='pow_screen']
  }else{
    if(!is.null(Data$fixed$pow_screen)){
      pow_screen = Data$fixed$pow_screen
    }else{
      pow_screen = -0.7
    }
  }
  if('alpha_SF' %in% names(parm)){
    parm[names(parm)=='alpha_SF']=.interval(parm[names(parm)=='alpha_SF'], -1, 0.6, reflect=FALSE)
    alpha_SF = 10^parm[names(parm)=='alpha_SF']
  }else{
    if(!is.null(Data$fixed$alpha_SF)){
      alpha_SF = Data$fixed$alpha_SF
    }else{
      alpha_SF = 1.5
    }
  }
  if('AGNfrac' %in% names(parm)){
    parm[names(parm)=='AGNfrac']=.interval(parm[names(parm)=='AGNfrac'], -2, -0.023, reflect=FALSE)
    AGNfrac = 10^parm[names(parm)=='AGNfrac']
  }else{
    if(!is.null(Data$fixed$AGNfrac)){
      AGNfrac = Data$fixed$AGNfrac
    }else{
      AGNfrac = 0
    }
  }
  
  # bad=1e10
  # if(any(as.numeric(massfunclist)[massfuncpos]<0)){if(verbose){print(-bad)}; return(bad)}
  # if(tau_birth < 0 | tau_birth>10){if(verbose){print(-bad)}; return(bad)}
  # if(tau_screen < 0 | tau_screen>10){if(verbose){print(-bad)}; return(bad)}
  # if(alpha_SF<0.1 | alpha_SF>4){if(verbose){print(-bad)}; return(bad)}
  # if(AGNfrac<0 | AGNfrac>0.95){if(verbose){print(-bad)}; return(bad)}
  
  # if(is.null(Dale)){
  #   Dale_Msol=NULL
  #   data('Dale_Msol', envir = environment())
  #   Dale=Dale_Msol
  # }
  
  # if(filters[1]=='all'){
  #   cenwave=NULL
  #   data('cenwave', envir = environment())
  #   filters=cenwave$filter
  #   filters=filters[filters %in% Data$flux$filter]
  #   Data$flux=Data$flux[Data$flux$filter %in% filters,]
  # }
  
  #SFH_dust=SFHfunc(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, z=redshift, outtype=NULL, speclib=Data$speclib)
  SFH_dust=do.call('SFHfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, tau_birth=tau_birth, tau_screen=tau_screen, pow_birth=pow_birth, pow_screen=pow_screen, z=redshift, outtype=NULL, speclib=Data$speclib, sparse=sparse), massfunclist))
  #if(sum(SFH_dust$massvec,na.rm=TRUE)==0){if(verbose){print(-bad)}; return(bad)}
  SFH_nodust=do.call('SFHfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, tau_birth=0, tau_screen=0, z=redshift, outtype=NULL, speclib=Data$speclib, sparse=sparse), massfunclist))
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale=Data$Dale)
  dustout=dustmass(Data$speclib$Wave[seq(1,dim(Data$speclib$Zspec[[1]])[2],by=sparse)], SFH_nodust$lum, SFH_dust$lum, Dale_in$Wave, Dale_in$Aspec)
  dustflux=Lum2Flux(Dale_in$Wave, Dale_in$Aspec*dustout[1], z=redshift)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  Dale_in=Dale_interp(alpha_SF=alpha_SF, AGNfrac=0, Dale=Data$Dale)
  Dale_in[,2]=Dale_in[,2]*dustout[3]
  Fractions=Dale_scale(alpha_SF=alpha_SF, AGNfrac=AGNfrac, Dale_in=Dale_in)
  #Might want to speed this next part up- spends much of the time loading the different filters!
  
  # if(is.null(filtout)){
  #   photom_out=photom_flux(finalspec[,1], finalspec[,2], outtype='Jy', filters = filters)$out
  # }else{
  #   photom_out={}
  #   fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  #   for(i in 1:length(filtout)){
  #     photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = filtout[[i]], lum = T)*1e23)
  #   }
  # }
  fluxnu=convert_wave2freq(finalspec[,2], finalspec[,1])
  photom_out={}
  for(i in 1:length(Data$filtout)){
    photom_out = c(photom_out, bandpass(flux = fluxnu, wave = finalspec[,1], filter = Data$filtout[[i]], lum = T)*1e23)
  }
  
  cutsig=(Data$flux$flux-photom_out)/Data$flux$fluxerr
  if(Data$like=='norm'){
    likeout=sum(dnorm(x=cutsig, log=TRUE), na.rm = TRUE)
  }else if(Data$like=='chisq'){
    likeout=dchisq(sum(cutsig^2), df=Data$N-length(parm), log=TRUE)
  }else if(Data$like=='st'){
    vardata = var(cutsig,na.rm = TRUE)
    dof=2*vardata/(vardata-1)
    #dof=interval(dof,0,Inf)
    dof=max(1, min(Inf, dof, na.rm = TRUE), na.rm = TRUE)
    likeout=sum(dt(cutsig, df=dof, log=TRUE), na.rm = TRUE)
    }else{
    stop('Bad like option!')
  }
  if(verbose){print(likeout)}
  if(Data$fit=='optim'){
    return(-likeout)
  }else if(Data$fit=='LD'){
    DustMass=log10(dustout[1])
    DustLum=log10(dustout[2])
    Dustfrac_bol=log10(Fractions[1])
    AGNfrac_bol=log10(Fractions[2])
    AGNLum=log10(dustout[2]*Fractions[2]/Fractions[1])
    Monitor=c(DustMass, DustLum, AGNLum, Dustfrac_bol, AGNfrac_bol)
    Monitor[!is.finite(Monitor)]=1e10
    return(list(LP=likeout,Dev=-2*likeout,Monitor=Monitor,yhat=1,parm=parm))
  }else if(Data$fit=='check'){
    SMtot=do.call('SMstarfunc', c(list(massfunc=massfunc, forcemass=forcemass, unimax=unimax, agescale=agescale, z=redshift), massfunclist))
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout, SMtot=SMtot))
  }
}

#This is a direct copy of the interval function from LaplacesDemon. Since I only use this one function I didn't want to

.interval=function (x, a = -Inf, b = Inf, reflect = TRUE) 
{
    if (missing(x)) 
        stop("The x argument is required.")
    if (a > b) 
        stop("a > b.")
    if (reflect & is.finite(a) & is.finite(b) & any(!is.finite(x))) {
        if (is.array(x)) {
            d <- dim(x)
            x <- as.vector(x)
        }
        x.inf.pos <- !is.finite(x)
        x[x.inf.pos] <- .interval(x[x.inf.pos], a, b, reflect = FALSE)
        if (is.array(x)) 
            x <- array(x, dim = d)
    }
    if (is.vector(x) & {
        length(x) == 1
    }) {
        if (reflect == FALSE) 
            x <- max(a, min(b, x))
        else if (x < a | x > b) {
            out <- TRUE
            while (out) {
                if (x < a) 
                  x <- a + a - x
                if (x > b) 
                  x <- b + b - x
                if (x >= a & x <= b) 
                  out <- FALSE
            }
        }
    }
    else if (is.vector(x) & {
        length(x) > 1
    }) {
        if (reflect == FALSE) {
            x.num <- which(x < a)
            x[x.num] <- a
            x.num <- which(x > b)
            x[x.num] <- b
        }
        else if (any(x < a) | any(x > b)) {
            out <- TRUE
            while (out) {
                x.num <- which(x < a)
                x[x.num] <- a + a - x[x.num]
                x.num <- which(x > b)
                x[x.num] <- b + b - x[x.num]
                if (all(x >= a) & all(x <= b)) 
                  out <- FALSE
            }
        }
    }
    else if (is.array(x)) {
        d <- dim(x)
        x <- as.vector(x)
        if (reflect == FALSE) {
            x.num <- which(x < a)
            x[x.num] <- a
            x.num <- which(x > b)
            x[x.num] <- b
        }
        else if (any(x < a) | any(x > b)) {
            out <- TRUE
            while (out) {
                x.num <- which(x < a)
                x[x.num] <- a + a - x[x.num]
                x.num <- which(x > b)
                x[x.num] <- b + b - x[x.num]
                if (all(x >= a) & all(x <= b)) 
                  out <- FALSE
            }
        }
        x <- array(x, dim = d)
    }
    return(x)
}
