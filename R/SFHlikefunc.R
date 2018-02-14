SFHlikefunc=function(parm=c(8,9,9,9,1,0.3,1.5,0), Data, massfit=c('burstmass', 'youngmass', 'oldmass', 'ancientmass'), taufit=c('tau_birth', 'tau_screen'), dustfit=c('alpha_SF', 'AGNfrac'), filters=c('FUV', 'NUV', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA', 'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1' , 'W2', 'W3', 'W4', 'P100', 'P160', 'S250' , 'S350', 'S500'), verbose=TRUE, like=TRUE){
  names(parm)=c(massfit, taufit, dustfit)
  if(verbose){print(parm)}
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
  if(burstmass<0){print(-bad); return(bad)}
  if(youngmass<0){print(-bad); return(bad)}
  if(oldmass<0){print(-bad); return(bad)}
  if(ancientmass<0){print(-bad); return(bad)}
  if(tau_birth < 0 | tau_birth>3){print(-bad); return(bad)}
  if(tau_screen < 0 | tau_screen>3){print(-bad); return(bad)}
  if(alpha_SF<0.1 | alpha_SF>4){print(-bad); return(bad)}
  if(AGNfrac<0 | AGNfrac>0.95){print(-bad); return(bad)}
  
  filters=filters[filters %in% Data$flux$filter]
  Data$flux=Data$flux[Data$flux$filter %in% filters,]
  
  SFH_dust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=tau_birth, tau_screen=tau_screen, z=Data$z, outtype='Jansky')
  SFH_nodust=SFHp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass, tau_birth=0, tau_screen=0, z=Data$z, outtype='Jansky')
  SMtot=SMstarp4(burstmass=burstmass, youngmass=youngmass, oldmass=oldmass, ancientmass=ancientmass)
  Dale_interp=Dale_interp(alpha_SF=alpha_SF, AGNfrac=AGNfrac)
  dustout=dustmass(SFH_dust$flux[,1], SFH_nodust$flux[,2], SFH_dust$flux[,2], Dale_Msol$Wave, Dale_interp, z=Data$z)
  dustflux=Lum2Flux(Dale_Msol$Wave, Dale_interp*dustout[1], z=Data$z)
  finalspec=addspec(dustflux[,1], dustflux[,2], SFH_dust$flux[,1], SFH_dust$flux[,2])
  photom_out=photom_flux(finalspec[,1], finalspec[,2], out='Jy', filters = filters)
  
  likeout=sum(dnorm(x=(Data$flux$flux-photom_out$out)/Data$flux$fluxerr, log=TRUE), na.rm = TRUE)
  if(verbose){print(likeout)}
  if(like){
    return(-likeout)
  }else{
    return(list(like=likeout, FinalSpec=finalspec, FinalPhotom=photom_out, SFH_dust=SFH_dust, SFH_nodust=SFH_nodust, Dust=dustout))
  }
}
