runShinySED <- function(flux = NULL, burstmass=1e8, youngmass=1e9, midmass=1e10, oldmass=1e10, ancientmass=1e10, 
                        z = 0.1, tau_birth = 1, tau_screen = 0.3, tau_AGN = 1,
                        alpha_SF_birth = 1, alpha_SF_screen = 3, alpha_SF_AGN = 0, AGNlum = 1e42) {
  
  shiny::shinyOptions(flux=flux, burstmass=log10(burstmass), youngmass=log10(youngmass),
               midmass=log10(midmass), oldmass=log10(oldmass), ancientmass=log10(ancientmass),
               z=z, tau_birth=tau_birth, tau_screen=tau_screen, tau_AGN=tau_AGN,
               alpha_SF_birth=alpha_SF_birth, alpha_SF_screen=alpha_SF_screen,
               alpha_SF_AGN=alpha_SF_AGN, AGNlum=log10(AGNlum))
  
  appDir <- system.file("shiny-examples", "app.R", package = "ProSpect")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  return(invisible(shiny::runApp(appDir, display.mode = "normal")))
}