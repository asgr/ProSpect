massfunc_const = function(age, mSFR=1, magemax=13.8e9){
  
  ifelse(age < magemax, mSFR, 0)

}

massfunc_p2=function(age, m1=1, m2=1, m1age=0, m2age=magemax, magemax=13.8e9){
  
  temp=rep(NA,length(age))
  temp[age<m1age]=0
  temp[age>m2age]=0
  sel=which(is.na(temp))
  temp[sel]= m1 + ((m2-m1)/(m2age-m1age))*(age[sel]-m1age)
  
  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}

massfunc_p3=function(age,m1=1, m2=1, m3=1, m1age=1e5, m2age=7e9, m3age=13e9, magemax=13.8e9){
  age[age<1e5]=1e5 #Since we make the spline in log space
  
  temp=splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age)

  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}

massfunc_p3_burst=function(age, mburst=0, m1=1, m2=1, m3=1, mburstage=1e8, m1age=1e5, m2age=7e9, m3age=13e9, magemax=13.8e9){
  age[age<1e5]=1e5 #Since we make the spline in log space
  
  temp=splinefun(c(m1age,m2age,m3age),c(m1,m2,m3),method='monoH.FC')(age)
  temp[age<mburstage]=temp[age<mburstage]+mburst
  
  temp[temp<0]=0
  temp[age>magemax]=0
  
  temp[age<mburstage]=temp[age<mburstage]+mburst
  
  invisible(temp)
}

massfunc_p4=function(age, m1=1, m2=1, m3=1, m4=1, m1age=1e5, m2age=2e9, m3age=9e9, m4age=13e9, magemax=13.8e9){
  age[age<1e5]=1e5 #Since we make the spline in log space
  
  temp=splinefun(log10(c(m1age,m2age,m3age,m4age)),c(m1,m2,m3,m4),method='monoH.FC')(log10(age))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}

massfunc_p6=function(age, m1=1, m2=1, m3=1, m4=1, m5=1, m6=1, m1age=1e5, m2age=1e8, m3age=1e9, m4age=5e9, m5age=9e9, m6age=1.3e10, magemax=13.8e9){
  age[age<1e5]=1e5 #Since we make the spline in log space
  
  temp=splinefun(log10(c(m1age,m2age,m3age,m4age,m5age,m6age)),c(m1,m2,m3,m4,m5,m6),method='monoH.FC')(log10(age))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}

massfunc_exp=function(age, mSFR=10, mtau=1, mpivot=magemax, magemax=13.8e9){
  
  temp=mSFR*exp(-mtau*((mpivot-age)/mpivot))
  
  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}

massfunc_CSFH=function(age, mSFR=10, mpeak=1e10, mperiod=1e9, mskew=0.5, magemax=13.8e9){
  
  temp=mSFR*dnorm(((age-mpeak)/mperiod)*(exp(mskew))^asinh((age-mpeak)/mperiod))*sqrt(2*pi)
  
  temp[temp<0]=0
  temp[age>magemax]=0
  invisible(temp)
}
