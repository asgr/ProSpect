speclib_download = function(stellpop = 'avail',
                            destpath = '',
                            type = 'rda',
                            URL = 'https://tinyurl.com/prospect-speclib/') {
  if (stellpop == 'avail') {
    url.show(paste0(URL, 'avail.txt?raw=1'))
  } else{
    speclib = paste0(stellpop, '.',tolower(type))
    
    download.file(paste0(URL, speclib, '?raw=1'),
                  destfile = paste0(destpath, speclib))
  }
}

speclib_FITSload = function(file, check = FALSE) {
  if (!requireNamespace("Rfits", quietly = TRUE)) {
    stop(
      'The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.',
      call. = FALSE
    )
  }
  
  output = list()
  temp = Rfits::Rfits_read(file,
                           header = FALSE,
                           pointer = FALSE,
                           data.table = FALSE)
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
  
  Zspec_loc = sort(grep(pattern = 'Zspec', names(temp), value = TRUE))
  Zspec = temp[Zspec_loc]
  names(Zspec) = NULL
  lapply(Zspec, as.matrix)
  output$Zspec = Zspec
  
  Zevo_loc = sort(grep(pattern = 'Zevo', names(temp), value = TRUE))
  Zevo = temp[Zevo_loc]
  names(Zevo) = NULL
  lapply(Zevo, as.data.frame)
  for (i in 1:length(Zevo)) {
    colnames(Zevo[[i]]) = c('SMstar', 'SMgas', 'SMtot', 'SFR', 'SMrem')
  }
  output$Zevo = Zevo
  
  if(check){
    speclib_check(output)
  }
  
  return(output)
}

speclib_check = function(speclib) {
  all_good = TRUE
  
  message(' - - - - ')
  message('General skeleton checks:')
  message(' - - - - ')
  
  temp = checkNumeric(
    speclib$Z,
    any.missing = FALSE,
    unique = TRUE,
    sorted = TRUE
  )
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Z check failed')
    print(temp)
  }else{
    message(' - Z check passed')
  }
  
  Z_len = length(speclib$Z)
  
  temp = checkNumeric(
    speclib$Age,
    any.missing = FALSE,
    unique = TRUE,
    sorted = TRUE
  )
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Age check failed')
    print(temp)
  }else{
    message(' - Age check passed')
  }
  
  Age_len = length(speclib$Age)
  
  temp = checkNumeric(speclib$AgeBins,
                      len = Age_len + 1,
                      any.missing = FALSE,
                      unique = TRUE,
                      sorted = TRUE)
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - AgeBins check failed')
    print(temp)
  }else{
    message(' - AgeBins check passed')
  }
  
  temp = checkNumeric(speclib$AgeWeights,
                      len = Age_len,
                      any.missing = FALSE)
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - AgeWeights check failed')
    print(temp)
  }else{
    message(' - AgeWeights check passed')
  }
  
  temp = checkNumeric(
    speclib$Wave,
    unique = TRUE,
    any.missing = FALSE,
    sorted = TRUE
  )
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Wave check failed')
    print(temp)
  }else{
    message(' - Wave check passed')
  }
  
  Wave_len = length(speclib$Wave)
  
  temp = checkList(speclib$Labels,
                   len = 5,
                   names = 'unique')
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Labels check failed')
    print(temp)
  }else{
    message(' - Labels check passed')
  }
  
  checkList(speclib$Zspec,
            len = Z_len,
            names = 'unnamed')
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Zspec check failed')
    print(temp)
  }else{
    message(' - Zspec check passed')
  }
  
  temp = checkList(speclib$Zevo,
                   len = Z_len,
                   names = 'unnamed')
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - Zevo check failed')
    print(temp)
  }else{
    message(' - Zevo check passed')
  }
  
  message(' - - - - ')
  message('Detailed checks:')
  message(' - - - - ')
  
  #Check Zspec contents
  for (i in 1:Z_len) {
    temp = checkMatrix(
      speclib$Zspec[[i]],
      mode = 'numeric',
      nrows = Age_len,
      ncols = Wave_len,
      any.missing = FALSE
    )
    
    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - - Zspec[[',i,']] check failed')
      print(temp)
    }else{
      message(' - - Zspec[[',i,']] check passed')
    }
  }
  
  #Check Zevo contents
  for (i in 1:Z_len) {
    temp = checkDataFrame(
      speclib$Zevo[[i]],
      types = 'numeric',
      nrows = Age_len,
      ncols = 5,
      any.missing = FALSE,
      col.names = 'unique'
    )
    
    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - - Zevo[[',i,']] data.frame check failed')
      print(temp)
    }else{
      message(' - - Zevo[[',i,']] data.frame check passed')
    }
    
    temp = checkNames(
      colnames(speclib$Zevo[[i]]),
      identical.to = c('SMstar', 'SMgas', 'SMtot', 'SFR', 'SMrem')
    )
    
    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - - Zevo[[',i,']] colnames check failed')
      print(temp)
    }else{
      message(' - - Zevo[[',i,']] colnames check passed')
    }
  }
  
  #Check Labels contents
  temp = checkNames(
    unlist(speclib$Labels, use.names = FALSE),
    identical.to = c(
      "Metallicity",
      "Time since ZAM / Yrs",
      "Wavelength / Ang",
      "Lsun / Ang (for 1 Msun SF)",
      "Lsun / Ang (for 1 Msun/Yr SFR)"
    )
  )
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - - Labels string check failed')
    print(temp)
  }else{
    message(' - - Labels string check passed')
  }
  
  temp = checkNames(
    names(speclib$Labels),
    identical.to = c("Zlab",
                     "Agelab",
                     "Wavelab",
                     "Lumlab",
                     "LumAgelab")
  )
  
  if (!isTRUE(temp)) {
    all_good = FALSE
    message(' - - Labels list.names check failed')
    print(temp)
  }else{
    message(' - - Labels list.names check passed')
  }
  
  message(' - - - - - - - - ')
  
  if (all_good) {
    message('Congratulations - all speclib checks are passing!')
  } else{
    message('Oh dear - not all speclib checks are passing :-(')
  }
  
  return(invisible(all_good))
}
