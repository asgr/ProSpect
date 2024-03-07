speclib_download = function(stellpop = 'avail',  destpath='', URL = 'https://tinyurl.com/prospect-speclib/'){
  if(stellpop == 'avail'){
    url.show(paste0(URL,'avail.txt?raw=1'))
  }else{
    speclib = paste0(stellpop, '.rda')
    
    download.file(paste0(URL, speclib, '?raw=1'),
                         destfile = paste0(destpath, speclib)
                         )
  }
}

load_speclib = function(file){
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
  
  output = list()
  temp = Rfits::Rfits_read(file, header=FALSE, pointer=FALSE, data.table=FALSE)
  output$Z = temp$Z
  output$Age = temp$Age
  output$AgeBins = temp$AgeBins
  output$AgeWeights = temp$AgeWeights
  output$Wave = temp$Wave
  
  output$Labels$Zlab = "Metallicity"
  output$Labels$Agelab = "Time since ZAM / Yrs"
  output$Labels$Wavelab = "Wavelength / Ang"
  output$Labels$Lumlab = "Lsun / Ang (for 1 Msun SF)"
  output$Labels$LumAgelab = "Lsun / Ang (for 1 Msun/Yr SFR)"
    
  Zspec_loc = grep(pattern = 'Zspec', names(temp))
  Zspec = temp[Zspec_loc]
  names(Zspec) = NULL
  output$Zspec = Zspec
  
  Zevo_loc = grep(pattern = 'Zevo', names(temp))
  Zevo = temp[Zevo_loc]
  names(Zevo) = NULL
  for(i in 1:length(Zevo)){
    colnames(Zevo[[i]]) = c('SMstar', 'SMgas', 'SMtot', 'SFR', 'SMrem')
  }
  output$Zevo = Zevo
  
  return(output)
}

check_specib = function(speclib){
  
}
