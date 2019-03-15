SFHp4=function(burstmass=1e8, youngmass=1e9, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), oldage=c(1e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7,  filters='all', Z=c(5,5,5,5), z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', cossplit=c(9e9,1.3e10), dosplit=FALSE, sparse=5, unimax=13.8e9, ...){
  
  burstmass=.interval(burstmass,0,Inf,reflect=FALSE)
  youngmass=.interval(youngmass,0,Inf,reflect=FALSE)
  oldmass=.interval(oldmass,0,Inf,reflect=FALSE)
  ancientmass=.interval(ancientmass,0,Inf,reflect=FALSE)
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(c(mean(burstage),mean(youngage),mean(oldage),mean(ancientage))),Z_args))
    Z=interp_param(Z,speclib$Z)$ID_mode
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in unique(Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }
  
  if(!is.null(filters)){
    if(filters[1]=='all'){
      cenwave=NULL
      data('cenwave', envir = environment())
      filters=cenwave$filter
    }
  }

  if(unimax!=FALSE & z>0){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    speclib$AgeWeights[speclib$Age > unimax-TravelTime]=0
  }
  
  if(length(Z)==1){Z=rep(Z,4)}
  if(dosplit){
    Tsplit=cossplit[1]-TravelTime
    Tstart=cossplit[2]-TravelTime
    oldage[2]=Tsplit
    ancientage[1]=Tsplit
    ancientage[2]=Tstart
  }
  
  speclib_burst=speclib$Zspec[[Z[1]]]
  speclib_young=speclib$Zspec[[Z[2]]]
  speclib_old=speclib$Zspec[[Z[3]]]
  speclib_ancient=speclib$Zspec[[Z[4]]]
  
  if(burstmass>0){
    burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
    burstage[1]=speclib$Age[burstageloc[1]]
    burstage[2]=speclib$Age[burstageloc[2]]
    burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
    burstlum[is.nan(burstlum)]=0
  }else{
    burstlum=0
  }
  
  if(youngmass>0){
    youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
    youngage[1]=speclib$Age[youngageloc[1]]
    youngage[2]=speclib$Age[youngageloc[2]]
    younglum=colSums(rbind(speclib_young[youngageloc[1]:youngageloc[2],])*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])*youngmass/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
    younglum[is.nan(younglum)]=0
  }else{
    younglum=0
  }
  
  if(oldmass>0){
    oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
    oldage[1]=speclib$Age[oldageloc[1]]
    oldage[2]=speclib$Age[oldageloc[2]]
    oldlum=colSums(rbind(speclib_old[oldageloc[1]:oldageloc[2],])*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])*oldmass/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
    oldlum[is.nan(oldlum)]=0
  }else{
    oldlum=0
  }
  
  if(ancientmass>0){
    ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))
    ancientage[1]=speclib$Age[ancientageloc[1]]
    ancientage[2]=speclib$Age[ancientageloc[2]]
    ancientlum=colSums(rbind(speclib_ancient[ancientageloc[1]:ancientageloc[2],])*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])*ancientmass/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
    ancientlum[is.nan(ancientlum)]=0
  }else{
    ancientlum=0
  }
  
  lum_unatten=burstlum+younglum+oldlum+ancientlum
  lum=lum_unatten
  
  lumtot_unatten=sum(c(0,diff(speclib$Wave))*lum)
  
  if(tau_birth!=0 & burstmass>0){
    lum=lum-burstlum
    speclib_burst[1:birthcloud,]=t(t(speclib_burst[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
    lum=lum+burstlum
    lumtot_birth=lumtot_unatten-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_birth=0
  }
  
  if(tau_screen!=0){
    lum=lum*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen)
    lumtot_screen=(lumtot_unatten-lumtot_birth)-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_screen=0
  }
  
  lumtot_atten=sum(c(0,diff(speclib$Wave))*lum)
  
  masstot=burstmass+youngmass+oldmass+ancientmass
  
  if(z<0 | is.null(filters)){
    return(invisible(list(wave_lum=speclib$Wave, lum_atten=lum, lum_unatten=lum_unatten,lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, masstot=masstot))) # returns the minimal luminosity outputs
  }
  if(z>0){
    flux=Lum2Flux(wave = speclib$Wave, lum = lum, z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }else{
    flux=cbind(wave = speclib$Wave, flux = lum*3e-07)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }

  ages=rbind(c(burstage, diff(burstage), mean(burstage)), c(youngage, diff(youngage), mean(youngage)), c(oldage, diff(oldage), mean(oldage)), c(ancientage, diff(ancientage), mean(ancientage)))
  colnames(ages)=c('lo','hi','duration','mean')
  ages=as.data.frame(ages)
  masses=data.frame(Forming=c(burstmass, youngmass, oldmass, ancientmass), Formed=c(burstmass+youngmass+oldmass+ancientmass,youngmass+oldmass+ancientmass, oldmass+ancientmass, ancientmass))
  SFR=masses$Forming/ages$duration
  sSFR=SFR/masses$Formed

  return(invisible(list(flux=flux, out=out, wave_lum=speclib$Wave, lum_unatten=lum_unatten, lum_atten=lum, lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, M2L=masstot/lumtot_unatten, ages=ages, masstot=masstot, masses=masses, SFR=SFR, sSFR=sSFR)))
}

SMstarp4=function(burstmass=1e8, youngmass=1e9, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), oldage=c(1e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=c(5,5,5,5), z=0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit=c(9e9,1.3e10), dosplit=FALSE, unimax=13.8e9, ...){
  
  burstmass=.interval(burstmass,0,Inf,reflect=FALSE)
  youngmass=.interval(youngmass,0,Inf,reflect=FALSE)
  oldmass=.interval(oldmass,0,Inf,reflect=FALSE)
  ancientmass=.interval(ancientmass,0,Inf,reflect=FALSE)
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(c(mean(burstage),mean(youngage),mean(oldage),mean(ancientage))),Z_args))
    Z=interp_param(Z,speclib$Z)$ID_mode
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(unimax!=FALSE & z>0){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    speclib$AgeWeights[speclib$Age > unimax-TravelTime]=0
  }
  
  if(length(Z)==1){Z=rep(Z,4)}
  if(dosplit){
    Tsplit=cossplit[1]-TravelTime
    Tstart=cossplit[2]-TravelTime
    oldage[2]=Tsplit
    ancientage[1]=Tsplit
    ancientage[2]=Tstart
  }

  burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))

  burstform=burstmass*sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]],na.rm=TRUE)/diff(burstage)
  youngform=youngmass*sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]],na.rm=TRUE)/diff(youngage)
  oldform=oldmass*sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]],na.rm=TRUE)/diff(oldage)
  ancientform=ancientmass*sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]],na.rm=TRUE)/diff(ancientage)
  
  totform=sum(burstform,youngform,oldform,ancientform,na.rm=TRUE)
  
  burststar=burstmass*sum(speclib$Zevo[[Z[1]]][burstageloc[1]:burstageloc[2],'SMstar']*speclib$AgeWeights[burstageloc[1]:burstageloc[2]],na.rm=TRUE)/diff(burstage)
  youngstar=youngmass*sum(speclib$Zevo[[Z[2]]][youngageloc[1]:youngageloc[2],'SMstar']*speclib$AgeWeights[youngageloc[1]:youngageloc[2]],na.rm=TRUE)/diff(youngage)
  oldstar=oldmass*sum(speclib$Zevo[[Z[4]]][oldageloc[1]:oldageloc[2],'SMstar']*speclib$AgeWeights[oldageloc[1]:oldageloc[2]],na.rm=TRUE)/diff(oldage)
  ancientstar=ancientmass*sum(speclib$Zevo[[Z[5]]][ancientageloc[1]:ancientageloc[2],'SMstar']*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]],na.rm=TRUE)/diff(ancientage)
  
  if(is.na(burststar) | is.nan(burststar)){burststar=0}
  if(is.na(youngstar) | is.nan(youngstar)){youngstar=0}
  if(is.na(oldstar) | is.nan(oldstar)){oldstar=0}
  if(is.na(ancientstar) | is.nan(ancientstar)){ancientstar=0}
  totstar=burststar+youngstar+oldstar+ancientstar
  return(c(BurstSMform=burstform, YoungSMform=youngform, OldSMform=oldform, AncientSMform=ancientform, BurstSMstar=burststar, YoungSMstar=youngstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMform=totform, TotSMstar=totstar))
}

SFHp5=function(burstmass=1e8, youngmass=1e9, midmass=1e10, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7, filters='all', Z=c(5,5,5,5,5), z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', cossplit=c(9e9,1.3e10), dosplit=FALSE, sparse=5, unimax=13.8e9, ...){
  
  burstmass=.interval(burstmass,0,Inf,reflect=FALSE)
  youngmass=.interval(youngmass,0,Inf,reflect=FALSE)
  midmass=.interval(midmass,0,Inf,reflect=FALSE)
  oldmass=.interval(oldmass,0,Inf,reflect=FALSE)
  ancientmass=.interval(ancientmass,0,Inf,reflect=FALSE)
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(c(mean(burstage),mean(youngage),mean(midage),mean(oldage),mean(ancientage))),Z_args))
    Z=interp_param(Z,speclib$Z)$ID_mode
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in unique(Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }
  
  if(!is.null(filters)){
    if(filters[1]=='all'){
      cenwave=NULL
      data('cenwave', envir = environment())
      filters=cenwave$filter
    }
  }
  
  if(unimax!=FALSE & z>0){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    speclib$AgeWeights[speclib$Age > unimax-TravelTime]=0
  }
  
  if(length(Z)==1){Z=rep(Z,5)}
  if(dosplit){
    Tsplit=cossplit[1]-TravelTime
    Tstart=cossplit[2]-TravelTime
    oldage[2]=Tsplit
    ancientage[1]=Tsplit
    ancientage[2]=Tstart
  }
  speclib_burst=speclib$Zspec[[Z[1]]]
  speclib_young=speclib$Zspec[[Z[2]]]
  speclib_mid=speclib$Zspec[[Z[3]]]
  speclib_old=speclib$Zspec[[Z[4]]]
  speclib_ancient=speclib$Zspec[[Z[5]]]
  
  if(burstmass>0){
    burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
    burstage[1]=speclib$Age[burstageloc[1]]
    burstage[2]=speclib$Age[burstageloc[2]]
    burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
    burstlum[is.nan(burstlum)]=0
  }else{
    burstlum=0
  }
  
  if(youngmass>0){
    youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
    youngage[1]=speclib$Age[youngageloc[1]]
    youngage[2]=speclib$Age[youngageloc[2]]
    younglum=colSums(rbind(speclib_young[youngageloc[1]:youngageloc[2],])*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])*youngmass/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
    younglum[is.nan(younglum)]=0
  }else{
    younglum=0
  }
  
  if(midmass>0){
    midageloc=c(which.min(abs(speclib$Age-midage[1])),which.min(abs(speclib$Age-midage[2])))
    midage[1]=speclib$Age[midageloc[1]]
    midage[2]=speclib$Age[midageloc[2]]
    midlum=colSums(rbind(speclib_mid[midageloc[1]:midageloc[2],])*speclib$AgeWeights[midageloc[1]:midageloc[2]])*midmass/sum(speclib$AgeWeights[midageloc[1]:midageloc[2]])
    midlum[is.nan(midlum)]=0
  }else{
    midlum=0
  }
  
  if(oldmass>0){
    oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
    oldage[1]=speclib$Age[oldageloc[1]]
    oldage[2]=speclib$Age[oldageloc[2]]
    oldlum=colSums(rbind(speclib_old[oldageloc[1]:oldageloc[2],])*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])*oldmass/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
    oldlum[is.nan(oldlum)]=0
  }else{
    oldlum=0
  }
  
  if(ancientmass>0){
    ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))
    ancientage[1]=speclib$Age[ancientageloc[1]]
    ancientage[2]=speclib$Age[ancientageloc[2]]
    ancientlum=colSums(rbind(speclib_ancient[ancientageloc[1]:ancientageloc[2],])*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])*ancientmass/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
    ancientlum[is.nan(ancientlum)]=0
  }else{
    ancientlum=0
  }
  
  lum=burstlum+younglum+midlum+oldlum+ancientlum
  lum_unatten=lum
  
  lumtot_unatten=sum(c(0,diff(speclib$Wave))*lum)
  
  if(tau_birth!=0 & burstmass>0){
    lum=lum-burstlum
    speclib_burst[1:birthcloud,]=t(t(speclib_burst[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
    lum=lum+burstlum
    lumtot_birth=lumtot_unatten-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_birth=0
  }
  
  if(tau_screen!=0){
    lum=lum*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen)
    lumtot_screen=(lumtot_unatten-lumtot_birth)-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_screen=0
  }
  
  lumtot_atten=sum(c(0,diff(speclib$Wave))*lum)
  
  masstot=burstmass+youngmass+midmass+oldmass+ancientmass
  
  if(z<0 | is.null(filters)){
    return(invisible(list(wave_lum=speclib$Wave, lum_atten=lum, lum_unatten=lum_unatten,lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, masstot=masstot))) # returns the minimal luminosity outputs
  }
  if(z>0){
    flux=Lum2Flux(wave = speclib$Wave, lum = lum, z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }else{
    flux=cbind(wave = speclib$Wave, flux = lum*3e-07)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }

  ages=rbind(c(burstage, diff(burstage), mean(burstage)), c(youngage, diff(youngage), mean(youngage)), c(midage, diff(midage), mean(midage)), c(oldage, diff(oldage), mean(oldage)), c(ancientage, diff(ancientage), mean(ancientage)))
  colnames(ages)=c('lo','hi','duration','mean')
  ages=as.data.frame(ages)
  masses=data.frame(Forming=c(burstmass, youngmass, midmass, oldmass, ancientmass), Formed=c(burstmass+youngmass+midmass+oldmass+ancientmass, youngmass+midmass+oldmass+ancientmass, midmass+oldmass+ancientmass, oldmass+ancientmass, ancientmass))
  SFR=masses$Forming/ages$duration
  sSFR=SFR/masses$Formed

  return(invisible(list(flux=flux, out=out, wave_lum=speclib$Wave, lum_unatten=lum_unatten, lum_atten=lum, lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, M2L=masstot/lumtot_unatten, ages=ages, masstot=masstot, masses=masses, SFR=SFR, sSFR=sSFR)))
}

SMstarp5=function(burstmass=1e8, youngmass=1e9, midmass=1e10, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=c(5,5,5,5,5), z=0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit=c(9e9,1.3e10), dosplit=FALSE, unimax=13.8e9, ...){
  
  burstmass=.interval(burstmass,0,Inf,reflect=FALSE)
  youngmass=.interval(youngmass,0,Inf,reflect=FALSE)
  midmass=.interval(midmass,0,Inf,reflect=FALSE)
  oldmass=.interval(oldmass,0,Inf,reflect=FALSE)
  ancientmass=.interval(ancientmass,0,Inf,reflect=FALSE)
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(c(mean(burstage),mean(youngage),mean(midage),mean(oldage),mean(ancientage))),Z_args))
    Z=interp_param(Z,speclib$Z)$ID_mode
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(unimax!=FALSE & z>0){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    speclib$AgeWeights[speclib$Age > unimax-TravelTime]=0
  }
  
  if(length(Z)==1){Z=rep(Z,5)}
  if(dosplit){
    Tsplit=cossplit[1]-TravelTime
    Tstart=cossplit[2]-TravelTime
    oldage[2]=Tsplit
    ancientage[1]=Tsplit
    ancientage[2]=Tstart
  }

  burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  midageloc=c(which.min(abs(speclib$Age-midage[1])),which.min(abs(speclib$Age-midage[2])))
  oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))

  burstform=burstmass*sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]],na.rm=TRUE)/diff(burstage)
  youngform=youngmass*sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]],na.rm=TRUE)/diff(youngage)
  midform=midmass*sum(speclib$AgeWeights[midageloc[1]:midageloc[2]],na.rm=TRUE)/diff(midage)
  oldform=oldmass*sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]],na.rm=TRUE)/diff(oldage)
  ancientform=ancientmass*sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]],na.rm=TRUE)/diff(ancientage)
  
  totform=sum(burstform,youngform,midform,oldform,ancientform,na.rm=TRUE)
  
  burststar=burstmass*sum(speclib$Zevo[[Z[1]]][burstageloc[1]:burstageloc[2],'SMstar']*speclib$AgeWeights[burstageloc[1]:burstageloc[2]],na.rm=TRUE)/diff(burstage)
  youngstar=youngmass*sum(speclib$Zevo[[Z[2]]][youngageloc[1]:youngageloc[2],'SMstar']*speclib$AgeWeights[youngageloc[1]:youngageloc[2]],na.rm=TRUE)/diff(youngage)
  midstar=midmass*sum(speclib$Zevo[[Z[3]]][midageloc[1]:midageloc[2],'SMstar']*speclib$AgeWeights[midageloc[1]:midageloc[2]],na.rm=TRUE)/diff(midage)
  oldstar=oldmass*sum(speclib$Zevo[[Z[4]]][oldageloc[1]:oldageloc[2],'SMstar']*speclib$AgeWeights[oldageloc[1]:oldageloc[2]],na.rm=TRUE)/diff(oldage)
  ancientstar=ancientmass*sum(speclib$Zevo[[Z[5]]][ancientageloc[1]:ancientageloc[2],'SMstar']*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]],na.rm=TRUE)/diff(ancientage)
  
  if(is.na(burststar) | is.nan(burststar)){burststar=0}
  if(is.na(youngstar) | is.nan(youngstar)){youngstar=0}
  if(is.na(midstar) | is.nan(midstar)){midstar=0}
  if(is.na(oldstar) | is.nan(oldstar)){oldstar=0}
  if(is.na(ancientstar) | is.nan(ancientstar)){ancientstar=0}
  
  totstar=burststar+youngstar+midstar+oldstar+ancientstar
  return(c(BurstSMform=burstform, YoungSMform=youngform, MidSMform=midform, OldSMform=oldform, AncientSMform=ancientform, BurstSMstar=burststar, YoungSMstar=youngstar, MidSMstar=midstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMform=totform, TotSMstar=totstar))
}

SFHfunc=function(massfunc=function(age, SFR=1){ifelse(age<1.3e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7, filters='all', Z=5, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', sparse=5, intSFR=FALSE, ...){
  
  dots=list(...)
  massfunc_args=dots[names(dots) %in% names(formals(massfunc))]
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(!is.null(filters)){
    if(filters[1]=='all'){
      cenwave=NULL
      data('cenwave', envir = environment())
      filters=cenwave$filter
    }
  }
  
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(speclib$Age*agescale),Z_args))
    Zlist=interp_param(Z, speclib$Z, log=TRUE)
    Zwmat=matrix(0, length(speclib$Age), length(speclib$Z))
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_hi)]=Zlist$weight_hi
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_lo)]=Zlist$weight_lo
    Zuse=which(colSums(Zwmat)>0)
    Zdoweight=TRUE
  }else{
    Zuse=Z
    Zdoweight=FALSE
  }
  
  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in Zuse){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }
  
  if(intSFR){
    massvec={}
    for(i in 1:length(speclib$Age)){
      tempint=try(integrate(massfunc, lower = speclib$AgeBins[i]*agescale, upper=speclib$AgeBins[i+1]*agescale)$value, silent=TRUE)
      if(class(tempint)=="try-error"){
        massvec=c(massvec,0)
      }else{
        massvec=c(massvec,tempint)
      }
    }
    massvec=massvec
  }else{
    massvec=do.call('massfunc',c(list(speclib$Age*agescale), massfunc_args))*speclib$AgeWeights
  }
  
  if(unimax!=FALSE & z>=0){
    agemax=unimax-cosdistTravelTime(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    massvec[speclib$Age>agemax]=0
  }
  if(forcemass==FALSE){
    masstot=sum(massvec)
  }else{
    masstot=sum(massvec)
    massvec=massvec*forcemass/masstot
    masstot=forcemass
  }
  
  if(length(Zuse)>1){
    lum=rep(0,length(speclib$Wave))
    for(Zid in Zuse){
      lum=lum+colSums(speclib$Zspec[[Zid]]*massvec*Zwmat[,Zid])
    }
  }else{
    lum=colSums(speclib$Zspec[[Zuse]]*massvec)
  }
  

  lumtot_unatten=sum(c(0,diff(speclib$Wave))*lum)
  lum_unatten=lum
  
  if(tau_birth!=0){
    lum=rep(0,length(speclib$Wave))
    for(Zid in Zuse){
      if(tau_birth!=0){
        speclib$Zspec[[Zid]][1:birthcloud,]=t(t(speclib$Zspec[[Zid]][1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
      }
      if(Zdoweight){
        lum=lum+colSums(speclib$Zspec[[Zid]]*massvec*Zwmat[,Zid])
      }else{
        lum=colSums(speclib$Zspec[[Zid]]*massvec)
      }
    }
    lumtot_birth=lumtot_unatten-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_birth=0
  }
  
  if(tau_screen!=0){
    lum=lum*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen)
    lumtot_screen=(lumtot_unatten-lumtot_birth)-sum(c(0,diff(speclib$Wave))*lum)
  }else{
    lumtot_screen=0
  }
  
  lumtot_atten=sum(c(0,diff(speclib$Wave))*lum)
  
  if(z<0 | is.null(filters)){
    return(invisible(list(wave_lum=speclib$Wave, lum_atten=lum, lum_unatten=lum_unatten,lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, masstot=masstot))) # returns the minimal luminosity outputs
  }
  if(z>0){
    flux=Lum2Flux(wave = speclib$Wave, lum = lum, z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }else{
    flux=cbind(wave = speclib$Wave, flux = lum*3e-07)
    if(!is.null(outtype)){
      out=photom_flux(flux, outtype = outtype, filters = filters)
      if(is.list(filters)){
        cenout={}
        for(i in filters){
          cenout=c(cenout,cenwavefunc(i))
        }
        out=data.frame(filter=NA, cenwave=cenout, out=out)
      }else{
        out=data.frame(cenwave[match(filters, cenwave$filter),], out=out)
      }
    }else{
      out=NULL
    }
  }
  
  return(list(flux=flux, out=out, wave_lum=speclib$Wave, lum_unatten=lum_unatten, lum_atten=lum, lumtot_unatten=lumtot_unatten, lumtot_atten=lumtot_atten, lumtot_birth=lumtot_birth, lumtot_screen=lumtot_screen, SFR=massvec/speclib$AgeWeights, masstot=masstot, massvec=massvec, M2L=masstot/lumtot_unatten))
}

SMstarfunc=function(massfunc=function(age, SFR=1){ifelse(age<1.3e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=5, z=0.1, H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM, ref, ...){
  
  dots=list(...)
  massfunc_args=dots[names(dots) %in% names(formals(massfunc))]
  
  if(stellpop=='BC03lr'){
    if(is.null(speclib)){
      BC03lr=NULL
      data('BC03lr', envir = environment())
      speclib=BC03lr
    }
  }
  if(stellpop=='BC03hr'){
    if(is.null(speclib)){
      BC03hr=NULL
      data('BC03hr', envir = environment())
      speclib=BC03hr
    }
  }
  if(stellpop=='EMILES'){
    if(is.null(speclib)){
      EMILES=NULL
      data('EMILES', envir = environment())
      speclib=EMILES
    }
  }
  
  if(any(speclib$Age<1e7)){
    birthcloud=max(which(speclib$Age<=1e7))
  }else{
    birthcloud=1
  }
  
  if(is.function(Z)){
    dots=list(...)
    Z_args=dots[names(dots) %in% names(formals(Z))]
    Z=do.call('Z',c(list(speclib$Age*agescale),Z_args))
    Zlist=interp_param(Z, speclib$Z, log=TRUE)
    Zwmat=matrix(0, length(speclib$Age), length(speclib$Z))
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_hi)]=Zlist$weight_hi
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_lo)]=Zlist$weight_lo
    Zuse=which(colSums(Zwmat)>0)
    Zdoweight=TRUE
  }else{
    Zuse=Z
    Zdoweight=FALSE
  }
  
  massvec=do.call('massfunc',c(list(speclib$Age*agescale), massfunc_args))*speclib$AgeWeights
  
  if(unimax!=FALSE & z>=0){
    TravelTime=min(ancientage[2],unimax-cosdistTravelTime(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9)
  }
  massvec[speclib$Age>TravelTime]=0
  if(forcemass!=FALSE){
    masstot=sum(massvec)
    massvec=massvec*forcemass/masstot
  }
  
  totstar=rep(0,length(massvec))
  
  for(Zid in Zuse){
    if(Zdoweight){
      totstar=totstar+speclib$Zevo[[Zid]][,'SMstar']*massvec*Zwmat[,Zid]
    }else{
      totstar=speclib$Zevo[[Zid]][,'SMstar']*massvec
    }
  }
  
  burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  midageloc=c(which.min(abs(speclib$Age-midage[1])),which.min(abs(speclib$Age-midage[2])))
  oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))
  
  burstrescale=(burstage[2]-burstage[1])/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  youngrescale=(youngage[2]-youngage[1])/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  midrescale=(midage[2]-midage[1])/sum(speclib$AgeWeights[midageloc[1]:midageloc[2]])
  oldrescale=(oldage[2]-oldage[1])/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientrescale=(ancientage[2]-ancientage[1])/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
  
  burstform=sum(massvec[burstageloc[1]:burstageloc[2]])*burstrescale
  burststar=sum(totstar[burstageloc[1]:burstageloc[2]])*burstrescale
  youngform=sum(massvec[youngageloc[1]:youngageloc[2]])*youngrescale
  youngstar=sum(totstar[youngageloc[1]:youngageloc[2]])*youngrescale
  midform=sum(massvec[midageloc[1]:midageloc[2]])*midrescale
  midstar=sum(totstar[midageloc[1]:midageloc[2]])*midrescale
  oldform=sum(massvec[oldageloc[1]:oldageloc[2]])*oldrescale
  oldstar=sum(totstar[oldageloc[1]:oldageloc[2]])*oldrescale
  ancientform=sum(massvec[ancientageloc[1]:ancientageloc[2]])*ancientrescale
  ancientstar=sum(totstar[ancientageloc[1]:ancientageloc[2]])*ancientrescale
  
  
  return(c(BurstSMform=burstform, YoungSMform=youngform, MidSMform=midform, OldSMform=oldform, AncientSMform=ancientform, BurstSMstar=burststar, YoungSMstar=youngstar, MidSMstar=midstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMform=sum(massvec,na.rm = TRUE), TotSMstar=sum(totstar,na.rm = TRUE)))
}
