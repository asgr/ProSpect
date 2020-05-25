.thermal_frac_1p4Ghz = function(alpha_mid = -0.753){
  nu1_nu2 = 325/1400
  return((nu1_nu2^alpha_mid - nu1_nu2^-0.8) / (nu1_nu2^-0.1 - nu1_nu2^-0.8))
}

RadioThermal_nu = function(freq=1.4, FIR=10e-14, Te=1e4, th_frac=0.1){
  #freq in GHz, FIR in W/m^2, Te in K
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  th_1p4Ghz = 1.4e10 * FIR * (Te/1e4)^0.45 * (1.4)^(-0.1) # Jy at 1.4 Ghz
  
  th_flux = th_1p4Ghz * (freq/1.4)^-0.1
  sy_flux = (th_1p4Ghz * (1 - th_frac) / th_frac)  * (freq/1.4)^-0.8
  tot_flux = th_flux + sy_flux
  
  output=matrix(nrow=length(freq), ncol=3)
  output[,1] = th_flux
  output[,2] = sy_flux
  output[,3] = tot_flux
  colnames(output) = c('th_flux', 'sy_flux', 'tot_flux')
  return(invisible(output))
}

RadioThermal_lam = function(wave=1e10*.c_to_mps/(1e9*1.4), FIR=10e-14, Te=1e4, th_frac=0.1){
  #wave in Ang, FIR in W/m^2, Te in K
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  return(RadioThermal_nu(freq=.c_to_mps/(wave/1e10)/1e9, FIR=FIR, Te=Te, th_frac=th_frac))
}

Radio_add = function(wave, flux_freq, z=0, Te=1e4, th_frac=0.1, freq_samp=seq(-1,3,by=0.01), wavecut=3e7){
  #wave in Ang, flux_freq in Jy
  #z reflects the stretching of the input wave
  freq = .c_to_mps/(wave/(1+z)/1e10) # wave in Ang -> freq in Hz
  tempfun = approxfun(freq, flux_freq)
  #FIR in integral in W/m^2 between 42.5e-6m and 122.5e-6m (in the FIR)
  FIR = integrate(tempfun, 2.447285e+12, 7.05394e+12)$value*1e-26 
  wave_samp = (1e10*.c_to_mps/(1e9*10^freq_samp))*(1+z)
  selwave = wave/(1+z) > 1e6 #only subtract off FIR and longer
  radio_sub = RadioThermal_lam(wave=wave[selwave]/(1+z), FIR=FIR, Te=Te, th_frac=0.104)[,'tot_flux']
  flux_freq[selwave] = flux_freq[selwave] - radio_sub
  flux_freq[flux_freq<=0] = 1e-200
  radio_add = RadioThermal_lam(wave=wave_samp, FIR=FIR, Te=Te, th_frac=th_frac)
  output = addspec(wave, flux_freq,
                   wave_samp, radio_add[,'tot_flux'],
                   extrap=0)
  return(invisible(output))
}
