SKIRTOR_interp = function(lum = 1e+44, ta = 1, p = 1, q = 1, ct = 40, rm = 60,
                             an=30, SKIRTOR = NULL){
  
  if(is.null(SKIRTOR)){
    data('SKIRTOR', envir = environment()) 
  }
  
  
  tmix = interp_quick(ta, SKIRTOR$ta)
  pmix = interp_quick(p, SKIRTOR$p)
  qmix = interp_quick(q, SKIRTOR$q)
  oamix = interp_quick(ct, SKIRTOR$ct)
  rmmix = interp_quick(rm, SKIRTOR$rm)
  imix = interp_quick(an, SKIRTOR$an)
 
  slice = SKIRTOR$Aspec[c(tmix[1:2]), c(pmix[1:2]), c(qmix[1:2]),
                      c(oamix[1:2]), c(rmmix[1:2]), 1, c(imix[1:2]), ]
   
  weights = rep(1, 64)
  weights = weights * rep(tmix[3:4], each = 1, times = 32)
  weights = weights * rep(pmix[3:4], each = 2, times = 16)
  weights = weights * rep(qmix[3:4], each = 4, times = 8)
  weights = weights * rep(oamix[3:4], each = 8, times = 4)
  weights = weights * rep(rmmix[3:4], each = 16, times = 2)
  weights = weights * rep(imix[3:4], each = 32, times = 1)
  
  tempmat = matrix(as.numeric(slice), 64, 132)
  agn_spectrum = (colSums(tempmat * weights))/SKIRTOR$Wave
  
  
  agn_spectrum = (agn_spectrum/(3.839e33 * 1E11))*lum
  return(data.frame(wave=SKIRTOR$Wave, lum = agn_spectrum)) 

}
