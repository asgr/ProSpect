Zfunc_p2=function(age, Z1=0.02, Z2=Z1, Z1age=0, Z2age=Zagemax, Zagemax=13.8e9){
  
  temp=rep(NA,length(age))
  temp[age<Z1age]=Z1
  temp[age>Z2age]=1e-04
  sel=which(is.na(temp))
  temp[sel]= Z1 + ((Z2-Z1)/(Z2age-Z1age))*(age[sel]-Z1age)
  
  temp[temp<1e-04]=1e-04
  temp[age>Zagemax]=1e-04
  invisible(temp)
}