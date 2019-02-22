runShinySED <- function() {
  appDir <- system.file("shiny-examples", "app.R", package = "ProSpect")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}