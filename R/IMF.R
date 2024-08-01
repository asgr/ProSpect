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

IMF_Chabrier = function(mass, alpha = 2.3, a = 0.08, b = 0.69, masslow = 0.01, massmax = 150, massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25){
  norm = integrate(.Chabrier_shape, lower=masslow, upper=massmax, alpha=alpha, a=a, b=b, masslow=masslow, massmax=massmax, massmult=TRUE, rel.tol=rel.tol)$value

  if(is.numeric(mass)){
    return(massform*.Chabrier_shape(mass=mass, alpha, a=a, b=b, masslow=masslow, massmax=massmax, massmult=massmult)/norm)
  }else if(is.list(mass)){
    if(is.null(mass$lo) | is.null(mass$hi)){
      stop('lo and hi list elements must be present!')
    }
    if(length(mass$lo) != length(mass$hi)){
      stop('lo and hi must be the same length!')
    }
    if(any(mass$lo > mass$hi)){
      stop('some lo > hi')
    }
    i = NULL
    tempout = rep(0,length(mass$lo))
    for(i in 1:length(mass$lo)){
      tempout[i] = integrate(.Chabrier_shape, lower=mass$lo[i], upper=mass$hi[i], alpha, a=a, b=b, masslow=masslow, massmax=massmax, massmult=massmult, rel.tol=rel.tol)$value
    }
    return(massform*tempout/norm)
  }
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

IMF_Kroupa = function(mass, alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3, masslow = 0.01, mass1 = 0.08, mass2 = 0.5, massmax = 150, massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25){
  norm=integrate(.Kroupa_shape, lower=masslow, upper=massmax, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=TRUE, rel.tol=rel.tol)$value

  if(is.numeric(mass)){
    return(massform*.Kroupa_shape(mass=mass, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult)/norm)
  }else if(is.list(mass)){
    if(is.null(mass$lo) | is.null(mass$hi)){
      stop('lo and hi list elements must be present!')
    }
    if(length(mass$lo) != length(mass$hi)){
      stop('lo and hi must be the same length!')
    }
    if(any(mass$lo > mass$hi)){
      stop('some lo > hi')
    }
    i = NULL
    tempout = rep(0,length(mass$lo))
    for(i in 1:length(mass$lo)){
      tempout[i] = integrate(.Kroupa_shape, lower=mass$lo[i], upper=mass$hi[i], alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, masslow = masslow, mass1 = mass1, mass2=mass2, massmax=massmax, massmult=massmult, rel.tol=rel.tol)$value
    }
    return(massform*tempout/norm)
  }
}

binlims = function(input, log=FALSE){
  N_in = length(input)

  if(log){
    input = log10(input)
  }

  in_diff = diff(input)

  bins = input[1:(N_in - 1L)] + in_diff/2
  bins = c(input[1] - in_diff[1]/2, bins, input[N_in] + in_diff[N_in - 1L]/2)

  if(log){
    bins = 10^bins
  }

  bin_lo = bins[1:N_in]
  bin_hi = bins[2:(N_in + 1L)]

  return(list(lo = bin_lo, hi = bin_hi))
}

