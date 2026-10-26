.thermal_frac_1p4Ghz = function(alpha_mid = -0.753){
  #not used any more
  nu1_nu2 = 325/1400
  return((nu1_nu2^alpha_mid - nu1_nu2^-0.8) / (nu1_nu2^-0.1 - nu1_nu2^-0.8))
}

#Free-free fraction (at 1.4 GHz) of the Dale templates, i.e. a 1:9 thermal:non-thermal mix. This is a
#fixed property of the templates, NOT a free parameter. The FIR is calibrated against the free-free radio,
#so this is the ff_frac at which the new radio must be generated to reproduce the template's own radio;
#radiocont uses it for the subtraction, and the addition uses the user's ff_frac.
.ff_frac_cal = 0.1

.radioemission_nu = function(freq=1.4e9, FIR=10e-14, Te=1e4, ff_frac=.ff_frac_cal, ff_power=-0.1, sy_power=-0.8){
  #freq in Hz, FIR in W/m^2, Te in K
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  #The FIR is calibrated against the free-free 1.4 GHz flux, so ff_flux is anchored on the FIR and the
  #synchrotron term is tied to it by the observed mix. This keeps the thermal term on its physical
  #prediction, at the cost of the TOTAL 1.4 GHz flux scaling as 1/ff_frac.
  ff_1p4Ghz = 13536775000 * FIR * (Te/1e4)^0.45 # Jy at 1.4 Ghz, first part is 1.4e10*(1.4^-0.1) = 13536775000

  ff_flux = ff_1p4Ghz * (freq/1.4e9)^ff_power
  sy_flux = ff_1p4Ghz * (1 - ff_frac)/ff_frac * (freq/1.4e9)^sy_power
  tot_flux = ff_flux + sy_flux

  output=matrix(nrow=length(freq), ncol=3)
  output[,1] = ff_flux
  output[,2] = sy_flux
  output[,3] = tot_flux
  colnames(output) = c('ff_flux', 'sy_flux', 'tot_flux')
  return(output)
}

.radioemission_lam = function(wave=2.141e9, FIR=10e-14, Te=1e4,
                              ff_frac=.ff_frac_cal, ff_power=-0.1, sy_power=-0.8){
  #wave in Ang, FIR in W/m^2, Te in K
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  return(.radioemission_nu(freq=.c_to_mps/(wave/1e10), FIR=FIR, Te=Te, ff_frac=ff_frac, ff_power=ff_power, sy_power=sy_power))
}

radiocont = function(wave, flux, z=0, Te=1e4, ff_frac=0.1, ff_power=-0.1, sy_power=-0.8,
                     wavesamp=seq(6,9.4,by=0.1), flux_in='freq', flux_out=flux_in, subtractonly = FALSE){
  #z reflects the stretching of the input wave
  wave = wave/(1 + z)
  flux = flux*(1 + z)
  #wave in Ang, flux in Jy
  if(!is.vector(wave)){
    if(dim(wave)[2]==2){
      flux = wave[,2]
      wave = wave[,1]
    }
  }
  freq = .c_to_mps/(wave/1e10) # wave in Ang -> freq in Hz
  if(flux_in == 'wave'){
    flux = convert_wave2freq(flux, wave)
  }
  #magplot(wave, flux, log='xy',type='l', xlim=range(c(wave,10^wavesamp)), ylim=c(1e-4,1e-3))
  tempfun = approxfun(log10(freq), log10(flux), yleft=-200, yright=-200)
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  #FIR = integrate(tempfun, 2.447285e+12, 7.05394e+12)$value * 1e-26 OLD LESS ACCURATE FIR
  #FIR = 0.1 * (10^tempfun(9.146128)) / (1.4e10 * (1.4)^(-0.1))
  FIR = 7.387284e-12 * 10^tempfun(9.146128) #same as above, just pre-computed
  #The 1.4 GHz pivot (log10 nu = 9.146128) anchors the whole FIR -> radio calibration, so it has to
  #fall inside the input grid. Outside it, tempfun returns the -200 sentinel and every radio term
  #silently vanishes, which would otherwise look like a legitimate zero-radio galaxy.
  pivot = 9.146128
  if(pivot < min(log10(freq)) | pivot > max(log10(freq))){
    warning('The 1.4 GHz pivot (log10 frequency = 9.146128, i.e. wavelength 2.14e9 Ang / 21.4 cm) ',
            'falls outside the input wavelength range [', sprintf('%.3f', min(log10(wave))), ', ',
            sprintf('%.3f', max(log10(wave))), ']. The FIR luminosity cannot be estimated, so no ',
            'radio continuum emission is subtracted or added. Extend the input (and any waveout ',
            'argument used to generate it) to at least log10 wavelength 9.331.', call. = FALSE)
  }
  #The subtraction below runs over all wave > 1e6, so the re-added emission has to cover that same
  #range. Otherwise everything beyond max(wavesamp) is zeroed by addspec(extrap=0), which produces a
  #cliff in the output spectrum. Extend the grid at its native spacing, then always terminate exactly
  #on max(wave) so the added emission reaches the end of the input grid.
  wavesamp_maxneed = max(log10(wave))
  if(wavesamp_maxneed > max(wavesamp)){
    wavesamp_dlog = diff(wavesamp)
    wavesamp_step = if(any(wavesamp_dlog > 0)) min(wavesamp_dlog[wavesamp_dlog > 0]) else 0.1
    wavesamp_nadd = floor((wavesamp_maxneed - max(wavesamp))/wavesamp_step)
    #Radio is smooth in log-log, so cap the tail resolution rather than risk a pathological blow-up
    wavesamp_nadd = min(wavesamp_nadd, 1e4)
    if(wavesamp_nadd >= 1){
      wavesamp = c(wavesamp, max(wavesamp) + wavesamp_step*seq_len(wavesamp_nadd))
    }
    wavesamp = 10^sort(unique(c(wavesamp, wavesamp_maxneed)))
  }else{
    wavesamp = 10^wavesamp
  }
  #freqsamp = .c_to_mps / (wavesamp / 1e10)
  selwave = wave > 1e6 #only subtract off FIR and longer
  #Both Te and ff_frac are left at their defaults (not passed down) because what is being removed is the
  #radio the Dale templates already contain, at their own 1e4 K and 1:9 thermal:non-thermal mix. The
  #addition below uses the user's Te and ff_frac. Getting the template's fraction right here matters: the
  #subtracted total is FF/ff_frac, so it must use .ff_frac_cal or the round trip will not be an identity.
  radio_sub = .radioemission_nu(freq=freq[selwave], FIR=FIR)[,'tot_flux']
  #lines(wave[selwave], radio_sub, col='blue',lwd=3)
  flux[selwave] = flux[selwave] - radio_sub
  flux[flux<=0] = 1e-200
  #Coerced rather than compared against the assignable F symbol; this also keeps numeric 0/1 working
  if(!isTRUE(as.logical(subtractonly))){
    ## adding new radio in
    radio_add = .radioemission_lam(wave=wavesamp, FIR=FIR, Te=Te, ff_frac=ff_frac, ff_power=ff_power, sy_power=sy_power)[,'tot_flux']
    #lines(wavesamp, radio_add, col='red',lwd=1)
    output = addspec(wave, flux,
                     wavesamp, radio_add,
                     extrap=0)
    #lines(output, col='green')
    #legend('topright', legend=c('in','out','rem','add'), col=c('black','green','blue','red'), lwd=c(1,1,3,1))
  } else {
    ## only using subtract mode
    output = data.frame(wave = wave, flux = flux)
  }
  output$wave = output$wave*(1 + z)
  output$flux = output$flux/(1 + z)
  if(flux_out == 'wave'){
    output = data.frame(wave=output$wave, flux=convert_freq2wave(output$flux, output$wave))
  }
  return(output)
}
