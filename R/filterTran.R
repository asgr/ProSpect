filterTranMags = function(mag_in, mag_out, return='all'){
  if(is.vector(mag_in)){
    Nbands=1
  }else{
    Nbands = dim(mag_in)[2]
    inputnames = colnames(mag_in)
    mag_in = as.matrix(mag_in)
  }
  if(!is.vector(mag_out)){'mag_out should be a vector!'}
  
  fits=list()
  ref_bands={}
  
  if(Nbands==1){
    tempcol = mag_out - mag_in
    tempmean = mean(tempcol)
    tempsd = sd(tempcol)
    fits = c(beta = tempmean, sigma = tempsd)
  }else if(Nbands>=2){
    for(i in 1:Nbands){
      for(j in 1:Nbands){
        if(! abs(i-j)==1){next}
        bluesel = min(i,j)
        redsel = max(i,j)
        tempcol = mag_in[,bluesel] - mag_in[,redsel]

        tempfit = lm(mag_out - mag_in[,i] ~ tempcol)$coefficients[2:1]
        sigma = sd(mag_in[,i] + tempfit[1]*tempcol + tempfit[2] - mag_out, na.rm = TRUE)
        tempfit=c(tempfit,sigma)
        names(tempfit) = c('alpha', 'beta', 'sigma')
        ref_bands = rbind(ref_bands, c(inputnames[i], inputnames[bluesel], inputnames[redsel]))
        fits=c(fits, list(tempfit))
        names(fits)[length(fits)] = paste(inputnames[i],' + alpha.(',inputnames[bluesel],' - ',inputnames[redsel],') + beta +/- sigma', sep='')
      }
    }
  }else{
    stop('Number of bands for mag_in must be >=1')
  }
  if(return=='all'){
    return(fits)
  }else if (return=='bestalpha'){
    bestvec=rep(0,length(fits))
    for(i in 1:length(fits)){
      bestvec[i] = abs(fits[[i]]['alpha'])
    }
  }else if (return=='bestbeta'){
    bestvec=rep(0,length(fits))
    for(i in 1:length(fits)){
      bestvec[i] = abs(fits[[i]]['beta'])
    }
  }else if (return=='bestsigma'){
    bestvec=rep(0,length(fits))
    for(i in 1:length(fits)){
      bestvec[i] = fits[[i]]['sigma']
    }
  }else if (return=='bestall'){
    bestvec=rep(0,length(fits))
    for(i in 1:length(fits)){
      bestvec[i] = 2*(fits[[i]]['alpha']^2) + fits[[i]]['beta']^2 + fits[[i]]['sigma']^2
    }
  }else{
    stop('return must be one of all, bestalpha, bestbeta, bestsigma or bestall!')
  }
  
  params=fits[which.min(bestvec)]
  ref_bands = ref_bands[which.min(bestvec),]
  
  return(
    invisible(
      list(
        params = params,
        ref_bands = ref_bands,
        predict = mag_in[,ref_bands[1]] + params[[1]][1]*(mag_in[,ref_bands[2]] - mag_in[,ref_bands[3]]) + params[[1]][2]
      )
    )
  )
}

filterTranBands = function(filt_in, filt_out, zrange=c(0,0.5), Nsamp=1e3, seed=666, return='all'){
  if(is.list(filt_in)){
    Nfilt_in = length(filt_in)
  }else{
    Nfilt_in = 1
  }
  if(is.null(names(filt_in)) & Nfilt_in>1){stop('filt_in must have filter names for sensical results!')}
  #if(is.null(names(filt_out))){stop('filt_out must have filter names for sensical results!')}
  BC03lr=NULL
  data('BC03lr', envir = environment())
  Dale_NormTot=NULL
  data('Dale_NormTot', envir = environment())
  set.seed(seed)
  params = cbind(runif(Nsamp), runif(Nsamp), runif(Nsamp), runif(Nsamp), runif(Nsamp), runif(Nsamp,zrange[1],zrange[2]))
  ranSFHs = matrix(0,Nsamp,Nfilt_in+1)
  if(Nfilt_in==1 & !is.list(filt_in)){
    filt_in=list(filt_in)
  }
  filtout = c(filt_in, list(filt_out))
  for(i in 1:Nsamp){
    ranSFHs[i,] = ProSpectSED(m1=params[i,1], m2=params[i,2], m3=params[i,3], m4=params[i,4], m5=params[i,5],
                    z=params[i,6], filtout=filtout, speclib=BC03lr, Dale=Dale_NormTot,
                    returnall = FALSE)
  }
  ranSFHs = data.frame(ranSFHs)
  if(length(names(filt_out))==1){
    colnames(ranSFHs) = c(names(filt_in), names(filt_out))
  }else{
    colnames(ranSFHs) = c(names(filt_in), names('Target'))
  }
  for(i in 1:dim(ranSFHs)[2]){
    ranSFHs[,i] = Jansky2magAB(ranSFHs[,i])
  }
  
  fits = filterTranMags(mag_in = ranSFHs[,1:(dim(ranSFHs)[2]-1)], mag_out = ranSFHs[,dim(ranSFHs)[2]], return=return)
  return(fits)
}
