CF_general=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_birth=function(wave, tau=1.0,pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}

CF_screen=function(wave, tau=0.3, pow=-0.7, pivot=5500){
  return=exp(-tau*(wave/pivot)^(pow))
}
