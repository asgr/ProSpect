Zfunc_p2=function(age, Z1=0.02, Z2=Z1, Z1age=0, Z2age=Zagemax, Zagemax=13.8, ...){
  #Scale functions ages to years
  Z1age = Z1age*1e9
  Z2age = Z2age*1e9
  Zagemax = Zagemax*1e9

  temp = rep(NA,length(age))
  temp[age < Z1age] = Z1
  temp[age > Z2age] = Z2
  sel = which(is.na(temp))
  temp[sel] = Z1 + ((Z2-Z1)/(Z2age-Z1age))*(age[sel]-Z1age)

  temp[temp < 1e-04] = 1e-04
  temp[age > Zagemax] = 1e-04
  temp[is.na(temp)] = 1e-04
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
  #Build the total with the same rectangle rule used for the CDF, so the CDF is
  #self-consistent and ends exactly at 0. Adaptive integrate() disagrees with
  #the rectangle sum by ~0.8% and can miss a narrow burst entirely, returning a
  #silently wrong normalisation that flattens the whole Z history.
  agegrid = seq(0,Zagemax+1e8,by=1e8)
  massvec = massfunc(agegrid, ...)
  masssum = sum(massvec*1e8)
  if(masssum > 0){
    massCDF = 1 - cumsum(massvec*1e8)/masssum
    massCDF[which(massCDF < 0)] = 0
    massCDF[which(massCDF > 1)] = 1
    tempfunc = approxfun(agegrid, massCDF)
    temp = Zstart + tempfunc(age)*(Zfinal - Zstart)
    temp[temp < Zstart] = Zstart
    temp[age > Zagemax] = Zstart
    temp[is.na(temp)] = Zstart
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
  #Same rectangle rule as the CDF below, see Zfunc_massmap_lin for why
  agegrid = seq(0,Zagemax+1e8,by=1e8)
  massvec = massfunc(agegrid, ...)
  masssum = sum(massvec*1e8)
  if(masssum > 0){
    massCDF = (1-cumsum(massvec*1e8)/masssum)*massfrac_final
    massCDF[which(massCDF < 0)] = 0
    massCDF[which(massCDF > 1)] = 1
    tempfunc = approxfun(agegrid,massCDF)
    temp= -yield*log(1-tempfunc(age))+Zstart
    temp[temp < Zstart] = Zstart
    temp[age > Zagemax] = Zstart
    temp[is.na(temp)] = Zstart
    return(temp)
  }else{
    return(rep(Zfinal, length(age)))
  }
}
