massfunc_const = function(age, mSFR=1, magemax=13.8){
  #Scale functions ages to years
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  return(ifelse(age < magemax, mSFR, 0))
}

massfunc_p2=function(age, m1=1, m2=m1, m1age=0, m2age=magemax, magemax=13.8){
  #Scale functions ages to years
  m1age=m1age*1e9
  m2age=m2age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=rep(NA,length(age))
  temp[age<m1age]=0
  temp[age>m2age]=0
  sel=which(is.na(temp))
  temp[sel]= m1 + ((m2-m1)/(m2age-m1age))*(age[sel]-m1age)
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_p3=function(age ,m1=1, m2=m1, m3=m2, m1age=1e-4, m2age=7, m3age=13, magemax=13.8){
  #Scale functions ages to years
  m1age=m1age*1e9
  m2age=m2age*1e9
  m3age=m3age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age)

  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_p3_burst=function(age, mburst=0, m1=1, m2=m1, m3=m2, m1age=1e-4, m2age=7, m3age=13, mburstage=0.1, magemax=13.8){
  #Scale functions ages to years
  mburstage=mburstage*1e9
  m1age=m1age*1e9
  m2age=m2age*1e9
  m3age=m3age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age)
  
  temp[temp<0]=0
  temp[age>magemax]=0
  
  temp[age<mburstage]=temp[age<mburstage]+mburst
  
  temp[temp<0]=0
  temp[age>magemax]=0
  
  return(temp)
}

massfunc_p4=function(age, m1=1, m2=m1, m3=m2, m4=m3, m1age=1e-4, m2age=2, m3age=9, m4age=13, magemax=13.8){
  #Scale functions ages to years
  m1age=m1age*1e9
  m2age=m2age*1e9
  m3age=m3age*1e9
  m4age=m4age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=splinefun(log10(c(m1age,m2age,m3age,m4age)),c(m1,m2,m3,m4),method='monoH.FC')(log10(age))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_p6=function(age, m1=1, m2=m1, m3=m2, m4=m3, m5=m4, m6=m5, m1age=1e-4, m2age=0.1, m3age=1, m4age=5, m5age=9, m6age=13, magemax=13.8){
  #Scale functions ages to years
  m1age=m1age*1e9
  m2age=m2age*1e9
  m3age=m3age*1e9
  m4age=m4age*1e9
  m5age=m5age*1e9
  m6age=m6age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=splinefun(log10(c(m1age,m2age,m3age,m4age,m5age,m6age)),c(m1,m2,m3,m4,m5,m6),method='monoH.FC')(log10(age))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_b5=function(age, m1=1, m2=m1, m3=m2, m4=m3, m5=m4, m1age=0, m2age=0.1, m3age=1, m4age=5, m5age=9, m6age=13, magemax=13.8){
  #Scale functions ages to years
  m1age=m1age*1e9
  m2age=m2age*1e9
  m3age=m3age*1e9
  m4age=m4age*1e9
  m5age=m5age*1e9
  m6age=m6age*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=
  ifelse(age>m1age & age<m2age, m1,
    ifelse(age>m2age & age<m3age, m2,
      ifelse(age>m2age & age<m3age, m2,
        ifelse(age>m3age & age<m4age, m3,
          ifelse(age>m4age & age<m5age, m4,
            ifelse(age>m5age & age<m6age, m5,
              0
            )
          )
        )
      )
    )
  )
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_exp=function(age, mSFR=10, mtau=1, mpivot=magemax, magemax=13.8){
  #Scale functions ages to years
  mpivot=mpivot*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=mSFR*exp(-mtau*((mpivot-age)/mpivot))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_exp_burst=function(age, mburst=0, mSFR=10, mtau=1, mpivot=magemax, mburstage=0.1, magemax=13.8){
  #Scale functions ages to years
  mpivot=mpivot*1e9
  mburstage=mburstage*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  temp=mSFR*exp(-mtau*((mpivot-age)/mpivot))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  
  temp[age<mburstage]=temp[age<mburstage]+mburst
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_snorm=function(age, mSFR=10, mpeak=10, mperiod=1, mskew=0.5, magemax=13.8){
  #Scale functions ages to years
  mpeak=mpeak*1e9
  mperiod=mperiod*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  #normal SFH
  Xtemp=(age-mpeak)/mperiod
  Ytemp=Xtemp*(exp(mskew)^asinh(Xtemp))
  temp=mSFR*exp(-(Ytemp^2)/2)
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_snorm_burst=function(age, mburst=0, mSFR=10, mpeak=10, mperiod=1, mskew=0.5, mburstage=0.1, magemax=13.8){
  #Scale functions ages to years
  mpeak=mpeak*1e9
  mperiod=mperiod*1e9
  mburstage=mburstage*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  #normal SFH
  Xtemp=(age-mpeak)/mperiod
  Ytemp=Xtemp*(exp(mskew)^asinh(Xtemp))
  temp=mSFR*exp(-(Ytemp^2)/2)
  #burst
  temp[age<mburstage]=temp[age<mburstage]+mburst
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_snorm_trunc=function(age, mSFR=10, mpeak=10, mperiod=1, mskew=0.5, magemax=13.8, mtrunc=2){
  #Scale functions ages to years
  mpeak=mpeak*1e9
  mperiod=mperiod*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  #normal SFH
  Xtemp=(age-mpeak)/mperiod
  Ytemp=Xtemp*(exp(mskew)^asinh(Xtemp))
  temp=mSFR*exp(-(Ytemp^2)/2)
  #truncation
  timetrunc = abs(magemax - mpeak)
  temp = temp * (1 - pnorm(age, mean = mpeak + timetrunc/mtrunc, sd = timetrunc/(2*mtrunc)))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_snorm_burst_trunc=function(age, mburst=0, mSFR=10, mpeak=10, mperiod=1, mskew=0.5, mburstage=0.1, magemax=13.8, mtrunc=2){
  #Scale functions ages to years
  mpeak=mpeak*1e9
  mperiod=mperiod*1e9
  mburstage=mburstage*1e9
  magemax=magemax*1e9
  
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  #normal SFH
  #normal SFH
  Xtemp=(age-mpeak)/mperiod
  Ytemp=Xtemp*(exp(mskew)^asinh(Xtemp))
  temp=mSFR*exp(-(Ytemp^2)/2)
  #burst
  temp[age<mburstage]=temp[age<mburstage]+mburst
  #truncation
  timetrunc = abs(magemax - mpeak)
  temp = temp * (1 - pnorm(age, mean = mpeak + timetrunc/mtrunc, sd = timetrunc/(2*mtrunc)))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}

massfunc_dtau=function(age, mSFR = 10,  mpeak = 10, mtau = 2, magemax = 13.8){
  mtau = mtau*1e9
  mpeak = mpeak*1e9
  magemax = magemax*1e9
  mtrunc = mtau / (magemax - mpeak)
  age[age<1e5]=1e5 #Stop dodgy very yound stellar pops forming
  
  x = rep(0, length(age))
  x[age<=mpeak] = (mtau + mpeak - age[age<=mpeak]) / mtau
  x[age>mpeak] = ((mtau/mtrunc) + mpeak - age[age>mpeak]) / (mtau/mtrunc)
  temp = x * exp(-x)
  temp = (mSFR/exp(-1)) * temp
  temp[temp<0]=0
  temp[age>magemax]=0
  return(temp)
}
