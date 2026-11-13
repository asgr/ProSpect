Fritz_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, al = 4, be = -0.5,
                        Fritz = NULL){
  if(is.null(Fritz)){
    Fritz=NULL
    data('Fritz', envir = environment())
  }

  ctmix = interp_quick(ct, Fritz$ct)
  almix = interp_quick(al, Fritz$al)
  bemix = interp_quick(be, Fritz$be)
  tamix = interp_quick(ta, Fritz$ta)
  rmmix = interp_quick(rm, Fritz$rm)
  anmix = interp_quick(an, Fritz$an)

  slice = Fritz$Aspec[
    c(ctmix[1:2]),
    c(almix[1:2]),
    c(bemix[1:2]),
    c(tamix[1:2]),
    c(rmmix[1:2]),
    c(anmix[1:2]),
  ]

  weights = rep(1,64)
  weights = weights * rep(ctmix[3:4], each=1, times=32)
  weights = weights * rep(almix[3:4], each=2, times=16)
  weights = weights * rep(bemix[3:4], each=4, times=8)
  weights = weights * rep(tamix[3:4], each=8, times=4)
  weights = weights * rep(rmmix[3:4], each=16, times=2)
  weights = weights * rep(anmix[3:4], each=32, times=1)
  tempmat = matrix(as.numeric(slice),64,178)
  return(data.frame(wave=Fritz$Wave, lum=lum * colSums(tempmat * weights)))
}

SKIRTOR_interp = function(lum = 1e+44, ct = 40, rm = 60, an = 30, ta = 1, p = 1, q = 1,
                          SKIRTOR = NULL){

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

Fritz_anav_interp = function(lum = 1e+44, ct = 40, rm = 60, ta = 1, al = 4, be = -0.5,
                             Fritz_anav = NULL){
  if(is.null(Fritz_anav)){
    Fritz_anav = NULL
    data('Fritz_anav', envir = environment())
    Fritz_anav = get('Fritz_anav', envir = environment())
  }
  if(!inherits(Fritz_anav, 'Fritz_anav')){
    stop('Fritz_anav must be a Fritz_anav object.', call. = FALSE)
  }

  ctmix = interp_quick(ct, Fritz_anav$ct)
  almix = interp_quick(al, Fritz_anav$al)
  bemix = interp_quick(be, Fritz_anav$be)
  tamix = interp_quick(ta, Fritz_anav$ta)
  rmmix = interp_quick(rm, Fritz_anav$rm)

  slice = Fritz_anav$Aspec[
    c(ctmix[1:2]),
    c(almix[1:2]),
    c(bemix[1:2]),
    c(tamix[1:2]),
    c(rmmix[1:2]),
  ]

  weights = rep(1,32)
  weights = weights * rep(ctmix[3:4], each=1, times=16)
  weights = weights * rep(almix[3:4], each=2, times=8)
  weights = weights * rep(bemix[3:4], each=4, times=4)
  weights = weights * rep(tamix[3:4], each=8, times=2)
  weights = weights * rep(rmmix[3:4], each=16, times=1)
  tempmat = matrix(as.numeric(slice), 32, length(Fritz_anav$Wave))
  return(data.frame(wave = Fritz_anav$Wave, lum = lum * colSums(tempmat * weights)))
}


SKIRTOR_anav_interp = function(lum = 1e+44, ct = 40, rm = 60, ta = 1, p = 1, q = 1,
                               SKIRTOR_anav = NULL){

  if(is.null(SKIRTOR_anav)){
    SKIRTOR_anav = NULL
    data('SKIRTOR_anav', envir = environment())
    SKIRTOR_anav = get('SKIRTOR_anav', envir = environment())
  }
  if(!inherits(SKIRTOR_anav, 'SKIRTOR_anav')){
    stop('SKIRTOR_anav must be a SKIRTOR_anav object.', call. = FALSE)
  }

  tmix = interp_quick(ta, SKIRTOR_anav$ta)
  pmix = interp_quick(p, SKIRTOR_anav$p)
  qmix = interp_quick(q, SKIRTOR_anav$q)
  oamix = interp_quick(ct, SKIRTOR_anav$ct)
  rmmix = interp_quick(rm, SKIRTOR_anav$rm)

  slice = SKIRTOR_anav$Aspec[
    c(tmix[1:2]),
    c(pmix[1:2]),
    c(qmix[1:2]),
    c(oamix[1:2]),
    c(rmmix[1:2]),
    1,
  ]

  weights = rep(1,32)
  weights = weights * rep(tmix[3:4], each=1, times=16)
  weights = weights * rep(pmix[3:4], each=2, times=8)
  weights = weights * rep(qmix[3:4], each=4, times=4)
  weights = weights * rep(oamix[3:4], each=8, times=2)
  weights = weights * rep(rmmix[3:4], each=16, times=1)
  tempmat = matrix(as.numeric(slice), 32, length(SKIRTOR_anav$Wave))
  agn_spectrum = colSums(tempmat * weights)/SKIRTOR_anav$Wave
  agn_spectrum = (agn_spectrum/(3.839e33 * 1E11))*lum
  return(data.frame(wave = SKIRTOR_anav$Wave, lum = agn_spectrum))

}
