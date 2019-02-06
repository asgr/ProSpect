SFHp4=function(burstmass=1e8, youngmass=1e9, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), oldage=c(1e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7,  filters='all', Z=c(5,5,5,5), z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', cossplit=c(9e9,1.3e10), dosplit=FALSE, sparse=1){
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
  
  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in unique(Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }
  
  if(filters[1]=='all'){
    cenwave=NULL
    data('cenwave', envir = environment())
    filters=cenwave$filter
  }
  
  if(length(Z)==1){Z=rep(Z,4)}
  if(dosplit){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
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
  if(tau_birth!=0){
    speclib_burst[1:birthcloud,]=t(t(speclib_burst[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_young[1:birthcloud,]=t(t(speclib_young[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_old[1:birthcloud,]=t(t(speclib_old[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_ancient[1:birthcloud,]=t(t(speclib_ancient[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
  }
  if(tau_screen!=0){
    speclib_burst=t(t(speclib_burst)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_young=t(t(speclib_young)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_old=t(t(speclib_old)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_ancient=t(t(speclib_ancient)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
  }
  burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))

  burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  younglum=colSums(rbind(speclib_young[youngageloc[1]:youngageloc[2],])*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])*youngmass/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  oldlum=colSums(rbind(speclib_old[oldageloc[1]:oldageloc[2],])*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])*oldmass/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientlum=colSums(rbind(speclib_ancient[ancientageloc[1]:ancientageloc[2],])*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])*ancientmass/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
  lum=burstlum+younglum+oldlum+ancientlum
  lumtot=sum(c(0,diff(speclib$Wave))*lum)
  
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

  masstot=burstmass+youngmass+oldmass+ancientmass

  ages=rbind(c(burstage, diff(burstage), mean(burstage)), c(youngage, diff(youngage), mean(youngage)), c(oldage, diff(oldage), mean(oldage)), c(ancientage, diff(ancientage), mean(ancientage)))
  colnames(ages)=c('lo','hi','duration','mean')
  masses=cbind(c(burstmass, youngmass, oldmass, ancientmass), c(burstmass+youngmass+oldmass+ancientmass,youngmass+oldmass+ancientmass, oldmass+ancientmass, ancientmass))
  colnames(masses)=c('Forming','Formed')
  SFR=masses[,'Forming']/ages[,'duration']
  sSFR=SFR/masses[,'Formed']

  return=list(flux=flux, lum=lum, out=out, masstot=masstot, lumtot=lumtot, M2L=masstot/lumtot, call=match.call(), ages=ages, masses=masses, SFR=SFR, sSFR=sSFR)
}

SMstarp4=function(burstmass=1e8, youngmass=1e9, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), oldage=c(1e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=c(5,5,5,5), z=0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit=c(9e9,1.3e10), dosplit=FALSE){
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
  
  if(length(Z)==1){Z=rep(Z,4)}
  if(dosplit){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
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

  burststar=burstmass*sum(speclib$Zevo[[Z[1]]][burstageloc[1]:burstageloc[2],'SMstar']*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  youngstar=youngmass*sum(speclib$Zevo[[Z[2]]][youngageloc[1]:youngageloc[2],'SMstar']*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  oldstar=oldmass*sum(speclib$Zevo[[Z[3]]][oldageloc[1]:oldageloc[2],'SMstar']*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientstar=ancientmass*sum(speclib$Zevo[[Z[4]]][ancientageloc[1]:ancientageloc[2],'SMstar']*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
  totstar=burststar+youngstar+oldstar+ancientstar
  return(c(BurstSMstar=burststar, YoungSMstar=youngstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMstar=totstar))
}

SFHp5=function(burstmass=1e8, youngmass=1e9, midmass=1e10, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7, filters='all', Z=c(5,5,5,5,5), z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', cossplit=c(9e9,1.3e10), dosplit=FALSE, sparse=1){
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
  
  if(sparse>1){
    sparse=seq(1,dim(speclib$Zspec[[1]])[2],by=sparse)
    for(i in unique(Z)){
      speclib$Zspec[[i]]=speclib$Zspec[[i]][,sparse]
    }
    speclib$Wave=speclib$Wave[sparse]
  }
  
  if(filters[1]=='all'){
    cenwave=NULL
    data('cenwave', envir = environment())
    filters=cenwave$filter
  }
  
  if(length(Z)==1){Z=rep(Z,5)}
  if(dosplit){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
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
  if(tau_birth!=0){
    speclib_burst[1:birthcloud,]=t(t(speclib_burst[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_young[1:birthcloud,]=t(t(speclib_young[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_mid[1:birthcloud,]=t(t(speclib_mid[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_old[1:birthcloud,]=t(t(speclib_old[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    # speclib_ancient[1:birthcloud,]=t(t(speclib_ancient[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
  }
  if(tau_screen!=0){
    speclib_burst=t(t(speclib_burst)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_young=t(t(speclib_young)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_mid=t(t(speclib_mid)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_old=t(t(speclib_old)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    speclib_ancient=t(t(speclib_ancient)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
  }
  burstageloc=c(which.min(abs(speclib$Age-burstage[1])),which.min(abs(speclib$Age-burstage[2])))
  youngageloc=c(which.min(abs(speclib$Age-youngage[1])),which.min(abs(speclib$Age-youngage[2])))
  midageloc=c(which.min(abs(speclib$Age-midage[1])),which.min(abs(speclib$Age-midage[2])))
  oldageloc=c(which.min(abs(speclib$Age-oldage[1])),which.min(abs(speclib$Age-oldage[2])))
  ancientageloc=c(which.min(abs(speclib$Age-ancientage[1])),which.min(abs(speclib$Age-ancientage[2])))

  burstlum=colSums(rbind(speclib_burst[burstageloc[1]:burstageloc[2],])*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])*burstmass/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  younglum=colSums(rbind(speclib_young[youngageloc[1]:youngageloc[2],])*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])*youngmass/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  midlum=colSums(rbind(speclib_mid[midageloc[1]:midageloc[2],])*speclib$AgeWeights[midageloc[1]:midageloc[2]])*midmass/sum(speclib$AgeWeights[midageloc[1]:midageloc[2]])
  oldlum=colSums(rbind(speclib_old[oldageloc[1]:oldageloc[2],])*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])*oldmass/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientlum=colSums(rbind(speclib_ancient[ancientageloc[1]:ancientageloc[2],])*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])*ancientmass/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
  lum=burstlum+younglum+midlum+oldlum+ancientlum
  lumtot=sum(c(0,diff(speclib$Wave))*lum)
  
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

  masstot=burstmass+youngmass+midmass+oldmass+ancientmass

  ages=rbind(c(burstage, diff(burstage), mean(burstage)), c(youngage, diff(youngage), mean(youngage)), c(midage, diff(midage), mean(midage)), c(oldage, diff(oldage), mean(oldage)), c(ancientage, diff(ancientage), mean(ancientage)))
  colnames(ages)=c('lo','hi','duration','mean')
  masses=cbind(c(burstmass, youngmass, midmass, oldmass, ancientmass), c(burstmass+youngmass+midmass+oldmass+ancientmass, youngmass+midmass+oldmass+ancientmass, midmass+oldmass+ancientmass, oldmass+ancientmass, ancientmass))
  colnames(masses)=c('Forming','Formed')
  SFR=masses[,'Forming']/ages[,'duration']
  sSFR=SFR/masses[,'Formed']

  return=list(flux=flux, lum=lum, out=out, masstot=masstot, lumtot=lumtot, M2L=masstot/lumtot, call=match.call(), ages=ages, masses=masses, SFR=SFR, sSFR=sSFR)
}

SMstarp5=function(burstmass=1e8, youngmass=1e9, midmass=1e10, oldmass=1e10, ancientmass=1e10, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=c(5,5,5,5,5), z=0, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, cossplit=c(9e9,1.3e10), dosplit=FALSE){
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
  
  if(length(Z)==1){Z=rep(Z,5)}
  if(dosplit){
    TravelTime=cosdistTravelTime(z=z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
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

  burststar=burstmass*sum(speclib$Zevo[[Z[1]]][burstageloc[1]:burstageloc[2],'SMstar']*speclib$AgeWeights[burstageloc[1]:burstageloc[2]])/sum(speclib$AgeWeights[burstageloc[1]:burstageloc[2]])
  youngstar=youngmass*sum(speclib$Zevo[[Z[2]]][youngageloc[1]:youngageloc[2],'SMstar']*speclib$AgeWeights[youngageloc[1]:youngageloc[2]])/sum(speclib$AgeWeights[youngageloc[1]:youngageloc[2]])
  midstar=midmass*sum(speclib$Zevo[[Z[3]]][midageloc[1]:midageloc[2],'SMstar']*speclib$AgeWeights[midageloc[1]:midageloc[2]])/sum(speclib$AgeWeights[midageloc[1]:midageloc[2]])
  oldstar=oldmass*sum(speclib$Zevo[[Z[4]]][oldageloc[1]:oldageloc[2],'SMstar']*speclib$AgeWeights[oldageloc[1]:oldageloc[2]])/sum(speclib$AgeWeights[oldageloc[1]:oldageloc[2]])
  ancientstar=ancientmass*sum(speclib$Zevo[[Z[5]]][ancientageloc[1]:ancientageloc[2],'SMstar']*speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])/sum(speclib$AgeWeights[ancientageloc[1]:ancientageloc[2]])
  totstar=burststar+youngstar+midstar+oldstar+ancientstar
  return(c(BurstSMstar=burststar, YoungSMstar=youngstar, MidSMstar=midstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMstar=totstar))
}

SFHfunc=function(massfunc=function(age, SFR=1){ifelse(age<1e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, stellpop='BC03lr', speclib=NULL, tau_birth=1.0, tau_screen=0.3, pow_birth=-0.7, pow_screen=-0.7, filters='all', Z=5, z = 0.1, H0 = 67.8, OmegaM = 0.308, OmegaL = 1 - OmegaM, ref, outtype='mag', sparse=1, intSFR=FALSE, ...){
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
  
  if(filters[1]=='all'){
    cenwave=NULL
    data('cenwave', envir = environment())
    filters=cenwave$filter
  }
  
  if(is.function(Z)){
    Zlist=interp_param(Z(speclib$Age*agescale), speclib$Z, log=TRUE)
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
    massvec=massvec/speclib$AgeWeights
  }else{
    massvec=massfunc(speclib$Age*agescale, ...)
  }
  
  if(unimax!=FALSE){
    agemax=unimax-cosdistTravelTime(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    massvec[speclib$Age>agemax]=0
  }
  if(forcemass==FALSE){
    masstot=sum(massvec*speclib$AgeWeights)
  }else{
    masstot=sum(massvec*speclib$AgeWeights)
    massvec=massvec*forcemass/masstot
    masstot=forcemass
  }
  
  lum=rep(0,length(speclib$Wave))
  
  for(Zid in Zuse){
    speclib_all=speclib$Zspec[[Zid]]
    
    if(tau_birth!=0){
      speclib_all[1:birthcloud,]=t(t(speclib_all[1:birthcloud,])*CF_birth(speclib$Wave, tau=tau_birth, pow=pow_birth))
    }
    if(tau_screen!=0){
      speclib_all=t(t(speclib_all)*CF_screen(speclib$Wave, tau=tau_screen, pow=pow_screen))
    }
    
    if(Zdoweight){
      lum=lum+colSums(speclib_all*speclib$AgeWeights*massvec*Zwmat[,Zid])
    }else{
      lum=colSums(speclib_all*speclib$AgeWeights*massvec)
    }
  }
  
  lumtot=sum(c(0,diff(speclib$Wave))*lum)
  
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
  
  return=list(flux=flux, lum=lum, out=out, massvec=massvec, masstot=masstot, lumtot=lumtot, M2L=masstot/lumtot, call=match.call())
}

SMstarfunc=function(massfunc=function(age, SFR=1){ifelse(age<1e+10,SFR,0)}, forcemass=FALSE, unimax=13.8e9, agescale=1, burstage=c(0,1e8), youngage=c(1e8,1e9), midage=c(1e9,5e9), oldage=c(5e9,9e9), ancientage=c(9e9,1.3e10), stellpop='BC03lr', speclib=NULL, Z=5, z=0.1, H0=67.8, OmegaM=0.308, OmegaL=1-OmegaM, ref, ...){
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
    Zlist=interp_param(Z(speclib$Age*agescale), speclib$Z, log=TRUE)
    Zwmat=matrix(0, length(speclib$Age), length(speclib$Z))
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_hi)]=Zlist$weight_hi
    Zwmat[cbind(1:length(speclib$Age),Zlist$ID_lo)]=Zlist$weight_lo
    Zuse=which(colSums(Zwmat)>0)
    Zdoweight=TRUE
  }else{
    Zuse=Z
    Zdoweight=FALSE
  }
  
  massvec=massfunc(speclib$Age*agescale, ...)*speclib$AgeWeights
  
  if(unimax!=FALSE){
    agemax=unimax-cosdistTravelTime(z = z, H0 = H0, OmegaM = OmegaM, OmegaL = OmegaL, ref = ref)*1e9
    massvec[speclib$Age>agemax]=0
  }
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
  
  burstform=sum(massvec[burstageloc])
  burststar=sum(totstar[burstageloc])
  youngform=sum(massvec[youngageloc])
  youngstar=sum(totstar[youngageloc])
  midform=sum(massvec[midageloc])
  midstar=sum(totstar[midageloc])
  oldform=sum(massvec[oldageloc])
  oldstar=sum(totstar[oldageloc])
  ancientform=sum(massvec[ancientageloc])
  ancientstar=sum(totstar[ancientageloc])
  
  
  return(c(BurstSMform=burstform, YoungSMform=youngform, MidSMform=midform, OldSMform=oldform, AncientSMform=ancientform, BurstSMstar=burststar, YoungSMstar=youngstar, MidSMstar=midstar, OldSMstar=oldstar, AncientSMstar=ancientstar, TotSMform=sum(massvec), TotSMstar=sum(totstar)))
}
