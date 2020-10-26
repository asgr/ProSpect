SpecLoad = function(filename,
                    AgeBins = NULL,
                    Labels = list(Zlab = 'Metallicity', 
                                  Agelab = 'Time since ZAM / Yrs', 
                                  Wavelab = 'Wavelength / Ang', 
                                  Lumlab = 'Lsun / Ang (for 1 Msun SF)',
                                  LumAgelab = 'Lsun / Ang (for 1 Msun/Yr SFR)'
                    )
){
  if(requireNamespace("Rfits", quietly = TRUE)){
    spec_FITS = Rfits::Rfits_read_all(filename, header=FALSE, data.table=FALSE)
    print(spec_FITS)
    
    Z = spec_FITS[['Z']][]
    Age = spec_FITS[['Age']][]
    Wave = spec_FITS[['Wave']][]
    
    if(is.null(AgeBins)){
      if(is.null(spec_FITS[['AgeBins']])){
        stop('Need AgeBins vector in FITS file, or function of vector passed into SpecLoad')
      }
      AgeBins = spec_FITS[['AgeBins']][]
    }else{
      if(is.function(AgeBins)){
        AgeBins = AgeBins(Age)
      }
    }
    
    if(length(Age) != (length(AgeBins) - 1)){
      stop('Age vector inconsistent with length of AgeBins: length(AgeBins) = length(Age) + 1')
    }
    
    Zspec_locs = grep('Zspec', names(spec_FITS))
    
    if(length(spec_FITS$Z) != length(Zspec_locs)){
      stop('Number of Metalicity libraries not consistent between Z and Zspec (there should be the same number)')
    }
    if(length(spec_FITS) != 5 + 2*length(Zspec_locs)){
      stop('Number of Metalicity libraries not consistent between Z/Zspec and Zevo  (there should be the same number)')
    }
    
    Zspec = list()
    for(i in Zspec_locs){
      Zspec = c(Zspec, list(spec_FITS[[i]][]))
      if(dim(Zspec[[i-5]])[1] != length(Age)){
        stop(paste0('Age vector inconsistent with dim(Zspec[[',i,']])[1]'))
      }
      if(dim(Zspec[[i-5]])[2] != length(Wave)){
        stop(paste0('Wave vector inconsistent with dim(Zspec[[',i,']])[2]'))
      }
    }
    
    Zevo_locs = grep('Zevo', names(spec_FITS))
    Zevo = list()
    for(i in Zevo_locs){
      Zevo = c(Zevo, list(spec_FITS[[i]][]))
      if(dim(Zevo[[i-5-length(Zspec_locs)]])[1] != length(spec_FITS[['Age']])){
        stop(paste0('Age vector inconsistent with dim(Zevo[[',i,']])[1]'))
      }
      if(dim(Zevo[[i-5-length(Zspec_locs)]])[2] != 5){
        stop(paste0('dim(Zevo[[',i,']])[2] should equal 5'))
      }
      if(!all(colnames(Zevo[[i-5-length(Zspec_locs)]]) == c('SMstar', 'SMgas', 'SMtot', 'SFR', 'SMrem'))){
        stop(paste0('Column names of Zevo[[',i,']] are not exactly: SMstar / SMgas / SMtot / SFR / SMrem'))
      }
    }
    
    spec_list = list(
      Z = Z,
      Age = Age,
      AgeBins = AgeBins,
      AgeWeights = diff(AgeBins),
      Wave = Wave,
      Labels = Labels,
      Zspec = Zspec,
      Zevo = Zevo
    )
    
    return(invisible(spec_list))
  }else{
    "Require Rfits to load target FITS spectral library!"
  }
}
