Zfunc_p2=function(age, Z1=0.02, Z2=Z1, Z1age=0, Z2age=Zagemax, Zagemax=13.8, ...){
  #Scale functions ages to years
  Z1age=Z1age*1e9
  Z2age=Z2age*1e9
  Zagemax=Zagemax*1e9
  
  temp=rep(NA,length(age))
  temp[age<Z1age]=Z1
  temp[age>Z2age]=Z2
  sel=which(is.na(temp))
  temp[sel]= Z1 + ((Z2-Z1)/(Z2age-Z1age))*(age[sel]-Z1age)
  
  temp[temp<1e-04]=1e-04
  temp[age>Zagemax]=1e-04
  temp[is.na(temp)]=1e-04
  return(temp)
}

Zfunc_massmap_lin=function(age, Zstart=1e-4, Zfinal=0.02, Zagemax=13.8, massfunc, ...){
  #Scale functions ages to years
  if(missing(Zagemax) & 'magemax' %in% names(list(...))){
    Zagemax = list(...)$magemax
  }
  Zagemax=Zagemax*1e9
  if(missing(massfunc)){
    stop('Need massfunc!')
  }
  masssum = integrate(massfunc, 0, Zagemax, ...)$value
  if(masssum > 0){
    massvec = massfunc(seq(0,Zagemax+1e8,by=1e8), ...)
    massCDF = 1 - cumsum(massvec*1e8)/masssum
    massCDF[which(massCDF < 0)] = 0
    massCDF[which(massCDF < 0)] = 1
    tempfunc = approxfun(seq(0,Zagemax+1e8,by=1e8), massCDF)
    temp = Zstart + tempfunc(age)*(Zfinal-Zstart)
    temp[temp<1e-04]=1e-04
    temp[age>Zagemax]=1e-04
    temp[is.na(temp)]=1e-04
    return(temp)
  }else{
    return(rep(Zfinal, length(age)))
  }
}

Zfunc_massmap_box=function(age, Zstart=1e-4, Zfinal=0.02, yield=0.03, Zagemax=13.8, massfunc, ...){
  #Scale functions ages to years
  massfrac_final=1-exp(-(Zfinal-Zstart)/yield)
  if(missing(Zagemax) & 'magemax' %in% names(list(...))){
    Zagemax = list(...)$magemax
  }
  Zagemax=Zagemax*1e9
  if(missing(massfunc)){
    stop('Need massfunc!')
  }
  masssum=integrate(massfunc, 0, Zagemax, ...)$value
  if(masssum > 0){
    massvec=massfunc(seq(0,Zagemax+1e8,by=1e8), ...)
    massCDF=(1-cumsum(massvec*1e8)/masssum)*massfrac_final
    massCDF[which(massCDF < 0)] = 0
    massCDF[which(massCDF < 0)] = 1
    tempfunc=approxfun(seq(0,Zagemax+1e8,by=1e8),massCDF)
    temp= -yield*log(1-tempfunc(age))+Zstart
    temp[temp<1e-04]=1e-04
    temp[age>Zagemax]=1e-04
    temp[is.na(temp)]=1e-04
    return(temp)
  }else{
    return(rep(Zfinal, length(age)))
  }
}
