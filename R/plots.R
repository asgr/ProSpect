plot.ProSpectSED = function(x,
                            xlim = c(1e3, 1e7),
                            ylim = 'auto',
                            xlab = 'Wavelength (Ang)',
                            ylab = 'auto',
                            grid = TRUE,
                            type = 'lum',
                            lwd_main = 5,
                            lwd_comp = 5,
                            ...) {
  if (type == 'lum') {
    if (ylim[1] == 'auto') {
      if(!is.null(x$Stars)){
        ylim = c(quantile(x$FinalLum[, 2], 0.45),
                 max(x$StarsUnAtten, na.rm = TRUE))
      }else{
        ylim = c(quantile(x$FinalLum[, 2], 0.45),
                 max(x$FinalLum[, 2], na.rm = TRUE))
      }
    }
    if (ylab[1] == 'auto') {
      ylab = 'Luminosity Density (Lsol/Ang)'
    }

    if(!is.null(x$Stars)){
      layout(rbind(1, 2), heights = c(0.7, 0.3))
      par(oma = c(3.1, 3.1, 1.1, 2.1))
      par(mar = c(0, 0, 0, 0))
    }else{
      layout(1)
      par(mar = c(3.1, 3.1, 1.1, 1.1))
    }

    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$FinalLum,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        ...
      )
    } else{
      plot(
        x$FinalLum,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        ...
      )
      points(x$Data, col='red', pch=16)
    }
    if(!is.null(x$Stars)){
      lines(x$StarsUnAtten,
            col = 'blue',
            lty = 2,
            lwd = lwd_comp)
      lines(x$StarsAtten, col = 'darkgreen', lwd = lwd_comp)
      lines(x$DustEmit, col = 'brown', lwd = lwd_comp)
    }

    if(!is.null(x$AGN)){
      lines(x$AGN, col = 'purple', lwd = lwd_comp)
    }

    legend(
      'topright',
      legend = c(
        'Total Lum',
        'Star Un-Atten',
        'Stars Atten',
        'Dust Emit',
        'AGN'
      ),
      col = c('black', 'blue', 'darkgreen', 'brown', 'purple'),
      lty = c(1, 2, 1, 1, 1),
      lwd = c(lwd_main, lwd_comp, lwd_comp, lwd_comp, lwd_comp)
    )

    if(is.null(x$Stars)){
      return(invisible(NULL))
    }

    par(mar = c(0, 0, 0, 0))
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        majorn = c(5,3)
      )
      par(usr = c(
        par()$usr[1:2],
        -max(x$Stars$Zvec, na.rm = TRUE) * 0.04,
        max(x$Stars$Zvec, na.rm = TRUE) * 1.04
      ))
      lines(x$Stars$agevec / 1e9,
            x$Stars$Zvec,
            col = 'red',
            lwd = 2)
      magicaxis::magaxis(4, col.axis = 'red', axis.col = 'red', majorn=3)
      legend(
        'bottomright',
        legend = c('SFR', 'Z'),
        col = c('black', 'red'),
        lty = 1,
        lwd = c(lwd_main, lwd_comp)
      )
    } else{
      plot(
        x$Stars$agevec / 1e9,
        x$Stars$SFR,
        xlab = 'Age (Gyr)',
        ylab = 'SFR (Msol/Yr)',
        type = 'l',
        lwd = lwd_main
      )
      par(usr = c(
        par()$usr[1:2],
        -max(x$Stars$Zvec, na.rm = TRUE) * 0.04,
        max(x$Stars$Zvec, na.rm = TRUE) * 1.04
      ))
      lines(x$Stars$agevec / 1e9,
            x$Stars$Zvec,
            col = 'red',
            lwd = 2)
      axis(4, col = 'red', col.axis = 'red')
      legend(
        'bottomright',
        legend = c('SFR', 'Z'),
        col = c('black', 'red'),
        lty = 1,
        lwd = c(lwd_main, lwd_comp)
      )
    }
  } else if (type == 'flux') {
    if (ylim[1] == 'auto') {
      ylim = quantile(x$FinalFlux[x$FinalFlux[, 'flux'] > 0, 'flux'], c(0.05, 1))
    }
    if (ylab[1] == 'auto') {
      ylab = 'Flux Density (Jansky)'
    }
    if (requireNamespace("magicaxis", quietly = TRUE)) {
      magicaxis::magplot(
        x$FinalFlux,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        grid = grid,
        ...
      )
    } else{
      plot(
        x$FinalFlux,
        log = 'xy',
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        type = 'l',
        lwd = lwd_main,
        ...
      )
    }
  } else{
    stop('type argument must be one of lum or flux!')
  }
}

plot.ProSpectSEDlike = function(x,
                                xlim = c(1e3, 1e7),
                                ylim = 'auto',
                                xlab = 'Wavelength (Ang)',
                                ylab = 'auto',
                                grid = TRUE,
                                type = 'flux',
                                ...) {
  if (type == 'flux') {
    if(isTRUE(x$Data$mode == 'spec') | is.null(x$Data$filtout)){
      data_mode = 'spec'
      photom = x$SEDout$Photom[,'flux']
      comp_type = 'l'
    }else if(isTRUE(x$Data$mode == 'photom') | isTRUE(x$Data$mode == 'both') | !is.null(x$Data$filtout)){
      data_mode = 'photom'
      photom = x$SEDout$Photom
      comp_type = 'p'
    }

    if(isTRUE(x$Data$mode == 'both')){
      data_mode = 'both'
    }

    if(data_mode == 'photom' | data_mode == 'spec'){
      layout(rbind(1, 2), heights = c(0.7, 0.3))
      par(oma = c(3.1, 3.1, 1.1, 1.1))
      par(mar = c(0, 0, 0, 0))
    }else{
      layout(1)
      par(mar = c(3.1, 3.1, 1.1, 1.1))
    }

    plot(
      x$SEDout,
      xlim = xlim,
      ylim = ylim,
      xlab = '',
      ylab = ylab,
      grid = grid,
      type = 'flux',
      ...
    )

    input_names = colnames(x$Data$flux)
    if('pivwave' %in% input_names){
      wavename = 'pivwave'
    }else if('cenwave' %in% input_names){
      wavename = 'cenwave'
    }else if('wave' %in% input_names){
      wavename = 'wave'
    }else{
      stop('wavelength column name not recognised, must')
    }

    if(data_mode == 'photom' | data_mode == 'both'){
      points(x$Data$flux[, c(wavename, 'flux')], pch = 16, col = 'red')
    }

    if(data_mode == 'both'){
      lines(x$Data$spec[,c('wave', 'flux')], col='grey', lwd=0.5)
    }

      if(data_mode == 'photom' | data_mode == 'both'){
        magicaxis::magerr(x$Data$flux[, wavename],
                          x$Data$flux[, 'flux'],
                          ylo = x$Data$flux[, 'fluxerr'],
                          col = 'red')
      }else if(data_mode == 'spec'){
        magicaxis::magerr(x$Data$flux[, wavename],
                          x$Data$flux[, 'flux'],
                          ylo = x$Data$flux[, 'fluxerr'],
                          col =  hsv(alpha=0.5),
                          poly = TRUE,
                          border = NA)
      }else if(data_mode == 'both'){
        magicaxis::magerr(x$Data$spec[, wavename],
                          x$Data$spec[, 'flux'],
                          ylo = x$Data$spec[, 'fluxerr'],
                          col =  hsv(alpha=0.5),
                          poly = TRUE,
                          border = NA)
      }

    legend('topleft', legend = paste('LP =', round(x$LP, 3)))

    if(data_mode == 'photom' | data_mode == 'spec'){
      resid = (x$Data$flux[, 'flux'] - photom) / x$Data$flux[, 'fluxerr']
      par(mar = c(0, 0, 0, 0))
      magicaxis::magplot(
        x$Data$flux[, wavename],
        (x$Data$flux[, 'flux'] - photom) / x$Data$flux[, 'fluxerr'],
        type = comp_type,
        pch = 16,
        col = 'red',
        grid = grid,
        log = 'x',
        xlim = xlim,
        ylim = c(-4, 4),
        xlab = xlab,
        ylab = '(Data - Model)/Error'
      )
    }

  } else if (type == 'lum') {
    plot(
      x$SEDout,
      xlim = xlim,
      ylim = ylim,
      xlab = '',
      ylab = ylab,
      grid = grid,
      type = 'lum',
      ...
    )
  }
}
