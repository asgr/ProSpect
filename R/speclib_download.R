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
