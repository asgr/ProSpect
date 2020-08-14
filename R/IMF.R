.Chabrier_shape=function(mass, alpha = 2.3, a = 0.08, b = 0.69, masslow = 0.01, massmax = 150, massmult = FALSE){
  #for singles a=0.08 and b=0.69
  #for binaries a=0.22 and b=0.57
  k=(1/(log(10)*1))*exp(-(log10(1)-log10(a))^2/(2*b^2))
  out=ifelse(mass<1, (1/(log(10)*mass))*exp(-(log10(mass)-log10(a))^2/(2*b^2)),
    k*(mass^-alpha)
  )
  out[mass<masslow | mass>massmax]=0
  if(massmult){out=out*mass}
  return(out)
}

IMF_Chabrier = function(mass, alpha = 2.3, a = 0.08, b = 0.69, masslow = 0.01, massmax = 150, massform = 1, massmult = FALSE){
  norm=integrate(.Chabrier_shape, lower=masslow, upper=massmax, alpha=alpha, a=a, b=b, masslow=masslow, massmax=massmax, massmult=TRUE)$value
  return(massform*.Chabrier_shape(mass=mass, alpha, a=a, b=b, masslow=masslow, massmax=massmax, massmult=massmult)/norm)
}

.Kroupa_shape=function(mass, alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3, masslow = 0.01, mass1 = 0.08, mass2=0.5, massmax=150, massmult = FALSE){
  k2=(mass1^-alpha1)/(mass1^-alpha2)
  k3=k2*(mass2^-alpha2)/(mass2^-alpha3)
  out=ifelse(mass<mass1, mass^-alpha1,
    ifelse(mass<mass2, k2*mass^-alpha2,
      k3*mass^-alpha3
    )
  )
  out[mass<masslow | mass>massmax]=0
  if(massmult){out=out*mass}
  return(out)
}

IMF_Kroupa = function(mass, alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3, masslow = 0.01, mass1 = 0.08, mass2 = 0.5, massmax = 150, massform = 1, massmult = FALSE){
  norm=integrate(.Kroupa_shape, lower=masslow, upper=massmax, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=TRUE)$value
  return(massform*.Kroupa_shape(mass=mass, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult)/norm)
}
