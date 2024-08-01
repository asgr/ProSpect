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
    return(paste0(destpath, speclib))
  }
}

speclib_FITSload = function(file, Labels = list(Zlab = 'Metallicity',
                                                Agelab = 'Time since ZAM / Yrs',
                                                Wavelab = 'Wavelength / Ang',
                                                Lumlab = 'Lsun / Ang (for 1 Msun SF)',
                                                LumAgelab = 'Lsun / Ang (for 1 Msun/Yr SFR)'),
                            check = FALSE) {
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
  output$Z = temp[['Z']]
  output$Age = temp[['Age']]
  output$AgeBins = temp[['AgeBins']]
  output$AgeWeights = temp[['AgeWeights']]
  output$Wave = temp[['Wave']]

  output$Labels = Labels

  Zspec_loc = grep(pattern = 'Zspec', names(temp), value = TRUE)
  Zspec = temp[Zspec_loc]
  names(Zspec) = NULL
  lapply(Zspec, as.matrix)
  output$Zspec = Zspec

  Zevo_loc = grep(pattern = 'Zevo', names(temp), value = TRUE)
  Zevo = temp[Zevo_loc]
  names(Zevo) = NULL
  lapply(Zevo, as.data.frame)
  for (i in 1:length(Zevo)) {
    colnames(Zevo[[i]]) = c('SMstar', 'SMgas', 'SMtot', 'SFR', 'SMrem')
  }
  output$Zevo = Zevo

  for(i in 1:length(output$Zspec)){
    loc = which(output$Zspec[[i]] < 0)
    output$Zspec[[i]][loc] = 0
  }

  if(check){
    speclib_check(output)
  }

  return(output)
}

speclib_check = function(speclib,
                         structure = TRUE,
                         coverage = TRUE,
                         Labels = list(Zlab = 'Metallicity',
                                                Agelab = 'Time since ZAM / Yrs',
                                                Wavelab = 'Wavelength / Ang',
                                                Lumlab = 'Lsun / Ang (for 1 Msun SF)',
                                                LumAgelab = 'Lsun / Ang (for 1 Msun/Yr SFR)')) {

  if(structure){

    all_good = TRUE

    message(' - - - - ')
    message('Skeleton checks:')
    message(' - - - - ')

    temp = checkList(
      speclib,
      any.missing = FALSE,
      unique = TRUE,
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' speclib list check failed')
      print(temp)
    }else{
      message(' speclib list check passed')
    }

    temp = checkNames(
      names(speclib),
      identical.to = c("Z", "Age", "AgeBins", "AgeWeights", "Wave", "Labels", "Zspec", "Zevo")
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' speclib list names check failed')
      print(temp)
    }else{
      message(' speclib list names check passed')
    }

    message(' - - - - ')
    message('General sub-list checks:')
    message(' - - - - ')

    temp = checkVector(
      speclib$Z,
      strict = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Z vector check failed')
      print(temp)
    }else{
      message(' - Z vector check passed')
    }

    temp = checkNumeric(
      speclib$Z,
      any.missing = FALSE,
      unique = TRUE,
      sorted = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Z numeric check failed')
      print(temp)
    }else{
      message(' - Z numeric check passed')
    }

    Z_len = length(speclib$Z)

    temp = checkVector(
      speclib$Age,
      strict = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Age vector check failed')
      print(temp)
    }else{
      message(' - Age vector check passed')
    }

    temp = checkNumeric(
      speclib$Age,
      any.missing = FALSE,
      unique = TRUE,
      sorted = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Age numeric check failed')
      print(temp)
    }else{
      message(' - Age numeric check passed')
    }

    Age_len = length(speclib$Age)

    temp = checkVector(
      speclib$AgeBins,
      strict = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - AgeBins vector check failed')
      print(temp)
    }else{
      message(' - AgeBins vector check passed')
    }

    temp = checkNumeric(speclib$AgeBins,
                        len = Age_len + 1,
                        any.missing = FALSE,
                        unique = TRUE,
                        sorted = TRUE)

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - AgeBins numeric check failed')
      print(temp)
    }else{
      message(' - AgeBins numeric check passed')
    }

    temp = checkVector(
      speclib$AgeWeights,
      strict = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - AgeWeights vector check failed')
      print(temp)
    }else{
      message(' - AgeWeights vector check passed')
    }

    temp = checkNumeric(speclib$AgeWeights,
                        len = Age_len,
                        any.missing = FALSE)

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - AgeWeights numeric check failed')
      print(temp)
    }else{
      message(' - AgeWeights numeric check passed')
    }

    temp = all.equal(diff(speclib$AgeBins), speclib$AgeWeights, tolerance=1e-3)

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - AgeWeights consistent with AgeBins check failed')
      print(temp)
    }else{
      message(' - AgeWeights consistent with AgeBins check passed')
    }

    temp = checkVector(
      speclib$Wave,
      strict = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Wave vector check failed')
      print(temp)
    }else{
      message(' - Wave vector check passed')
    }

    temp = checkNumeric(
      speclib$Wave,
      unique = TRUE,
      any.missing = FALSE,
      sorted = TRUE
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Wave numeric check failed')
      print(temp)
    }else{
      message(' - Wave numeric check passed')
    }

    Wave_len = length(speclib$Wave)

    temp = checkList(speclib$Labels,
                     len = 5,
                     names = 'unique')

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Labels list check failed')
      print(temp)
    }else{
      message(' - Labels list check passed')
    }

    checkList(speclib$Zspec,
              len = Z_len,
              names = 'unnamed')

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Zspec list check failed')
      print(temp)
    }else{
      message(' - Zspec list check passed')
    }

    temp = checkList(speclib$Zevo,
                     len = Z_len,
                     names = 'unnamed')

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - Zevo list check failed')
      print(temp)
    }else{
      message(' - Zevo list check passed')
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
        message(' - - Zspec[[',i,']] matrix check failed')
        print(temp)
      }else{
        message(' - - Zspec[[',i,']] matrix check passed')
      }
    }

    for (i in 1:Z_len) {
      temp = all(speclib$Zspec[[i]] >= 0)

      if (!isTRUE(temp)) {
        all_good = FALSE
        message(' - - Zspec[[',i,']] +ve flux check failed')
        print(temp)
      }else{
        message(' - - Zspec[[',i,']] +ve flux check passed')
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

      temp = checkNumeric(speclib$Zevo[[i]][,'SMstar'],
                          lower = 0,
                          upper = 1,
                          any.missing = FALSE
                          )

      if (!isTRUE(temp)) {
        all_good = FALSE
        message(' - - Zevo[[',i,']]$SMstar numeric check failed')
        print(temp)
      }else{
        message(' - - Zevo[[',i,']]$SMstar numeric check passed')
      }
    }

    #Check Labels contents
    temp = checkNames(
      unlist(speclib$Labels, use.names = FALSE),
      identical.to = unlist(Labels, use.names = FALSE)
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
      identical.to = names(Labels)
    )

    if (!isTRUE(temp)) {
      all_good = FALSE
      message(' - - Labels list names check failed')
      print(temp)
    }else{
      message(' - - Labels list names check passed')
    }

    message(' - - - - - - - - ')

    if (all_good) {
      message('Congratulations - all speclib checks are passing!')
    }else{
      message('Oh dear - not all speclib checks are passing :-(')
    }

    message(' - - - - - - - - \n')
  }else{
    all_good = NULL
  }

  if(coverage){
    message(' - - - - - - - - ')

    Zinfo = c(min(speclib$Z), median(speclib$Z), max(speclib$Z))
    message('Z min / med / max: ',Zinfo[1],' ',Zinfo[2],' ',Zinfo[3],' [',length(speclib$Z),']')

    Ageinfo = c(min(speclib$Age), median(speclib$Age), max(speclib$Age))
    message('Age min / med / max: ',Ageinfo[1],' ',Ageinfo[2],' ',Ageinfo[3],' [',length(speclib$Age),']')

    Waveinfo = c(min(speclib$Wave), median(speclib$Wave), max(speclib$Wave))
    message('Wave min / med / max: ',Waveinfo[1],' ',Waveinfo[2],' ',Waveinfo[3],' [',length(speclib$Wave),']')

    WaveBin = diff(speclib$Wave)
    WaveBininfo = c(min(WaveBin), median(WaveBin), max(WaveBin))
    message('Wave Bin min / med / max: ',WaveBininfo[1],' ',WaveBininfo[2],' ',WaveBininfo[3])

    message(' - - - - - - - - \n')
  }

  return(invisible(all_good))
}
