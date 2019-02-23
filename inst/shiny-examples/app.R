BC03lr=NULL
data('BC03lr', envir = environment())
Dale_NormTot=NULL
data('Dale_NormTot', envir = environment())
AGN_UnOb_Sparse=NULL
data('AGN_UnOb_Sparse', envir = environment())
cenwave=NULL
data('cenwave', envir = environment())
filters=c('FUV_GALEX', 'NUV_GALEX', 'u_SDSS', 'g_SDSS', 'r_SDSS', 'i_SDSS', 'Z_VISTA',
'Y_VISTA', 'J_VISTA', 'H_VISTA', 'K_VISTA', 'W1_WISE' , 'W2_WISE', 'W3_WISE', 'W4_WISE',
'P100_Herschel', 'P160_Herschel', 'S250_Herschel' , 'S350_Herschel', 'S500_Herschel')
shortnames=c('NUV','FUV')
cenwave=cenwave[cenwave$filter %in% filters,]
filtout={}
for(i in cenwave$filter){filtout=c(filtout,list(getfilt(i)))}

.checkoption=function(option,default){
  if(!is.null(getShinyOption(option))){
    return(getShinyOption(option))
  }else{
    return(default)
  }
}

flux=.checkoption('flux',NULL)
z=.checkoption('z',0.1)
burstmass=.checkoption('burstmass',8)
youngmass=.checkoption('youngmass',9)
midmass=.checkoption('midmass',10)
oldmass=.checkoption('oldmass',10)
ancientmass=.checkoption('ancientmass',10)
AGNlum=.checkoption('AGNlum',42)
tau_birth=.checkoption('taubirth',1)
tau_screen=.checkoption('tauscreen',0.3)
tau_AGN=.checkoption('tauscreen',1)
alpha_SF_birth=.checkoption('alpha_SF_birth',1)
alpha_SF_screen=.checkoption('alpha_SF_screen',3)
alpha_SF_AGN=.checkoption('alpha_SF_AGN',0)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("ProSpect SED"),
   
   # Sidebar with a slider input for number of bins 
        column(2,
          numericInput("z",
                       "Redshift",
                       value = z,
                       min = 0,
                       max = Inf
            
          ),
          sliderInput("burstmass",
                     "Burst Mass:",
                     min = 1,
                     max = 12,
                     value = burstmass, step=0.1),
         sliderInput("youngmass",
                     "Young Mass:",
                     min = 1,
                     max = 12,
                     value = youngmass, step=0.1),
         sliderInput("midmass",
                     "Mid Mass:",
                     min = 1,
                     max = 12,
                     value = midmass, step=0.1),
         sliderInput("oldmass",
                     "Old Mass:",
                     min = 1,
                     max = 12,
                     value = oldmass, step=0.1),
         sliderInput("ancientmass",
                     "Ancient Mass:",
                     min = 1,
                     max = 12,
                     value = ancientmass, step=0.1),
         sliderInput("AGNlum",
                     "AGN Luminosity:",
                     min = 40,
                     max = 47,
                     value = AGNlum, step=0.1)
        ),
        column(2,
          actionButton("stop",
                       "Stop!"
          ),
         sliderInput("tau_birth",
                     "Tau Birth:",
                     min = 0,
                     max = 2,
                     value = tau_birth, step=0.1),
         sliderInput("tau_screen",
                     "Tau Screen:",
                     min = 0,
                     max = 2,
                     value = tau_screen, step=0.1),
         sliderInput("tau_AGN",
                     "Tau AGN:",
                     min = 0,
                     max = 2,
                     value = tau_AGN, step=0.1),
         sliderInput("alpha_SF_birth",
                     "Alpha SF Birth:",
                     min = 0,
                     max = 4,
                     value = alpha_SF_birth, step=0.1),
         sliderInput("alpha_SF_screen",
                     "Alpha SF Screen:",
                     min = 0,
                     max = 4,
                     value = alpha_SF_screen, step=0.1),
         sliderInput("alpha_SF_AGN",
                     "Alpha SF AGN:",
                     min = 0,
                     max = 4,
                     value = alpha_SF_AGN, step=0.1)
        ),
      
      # Show a plot of the generated distribution
      column(8, 
        tabsetPanel(type='tabs',
          tabPanel('Flux',
            mainPanel(
              plotOutput("SED_flux_plot", height="600px")
            )
          ),
          tabPanel('Luminosity',
            mainPanel(
              plotOutput("SED_lum_plot", height="450px"),
              plotOutput("SFH_plot", height="250px")
            )
          )
        )
      )
)

server <- function(input, output) {
  
  observe({
    if(input$stop > 0){
      stopApp(returnValue = list(
                  z=input$z,
                  tau_birth=input$tau_birth,
                  tau_screen=input$tau_screen,
                  tau_AGN=input$tau_AGN,
                  alpha_SF_birth=input$alpha_SF_birth,
                  alpha_SF_screen=input$alpha_SF_screen,
                  alpha_SF_AGN=input$alpha_SF_AGN,
                  AGNlum=10^input$AGNlum,
                  burstmass=10^input$burstmass,
                  youngmass=10^input$youngmass,
                  midmass=10^input$midmass,
                  oldmass=10^input$oldmass,
                  ancientmass=10^input$ancientmass
                )
          )
    }
  })

   output$SED_flux_plot <- renderPlot({
     
     SED=ProSpectSED(SFH=SFHp5,
                  z=input$z,
                  tau_birth=input$tau_birth,
                  tau_screen=input$tau_screen,
                  tau_AGN=input$tau_AGN,
                  alpha_SF_birth=input$alpha_SF_birth,
                  alpha_SF_screen=input$alpha_SF_screen,
                  alpha_SF_AGN=input$alpha_SF_AGN,
                  AGNlum=10^input$AGNlum,
                  speclib=BC03lr,
                  Dale=Dale_NormTot,
                  AGN=AGN_UnOb_Sparse,
                  filtout=filtout,
                  burstmass=10^input$burstmass,
                  youngmass=10^input$youngmass,
                  midmass=10^input$midmass,
                  oldmass=10^input$oldmass,
                  ancientmass=10^input$ancientmass
                  )
     
      par(mar=c(3.1,3.1,1.1,1.1))
     
      magicaxis::magplot(SED$FinalFlux,
              log='xy',
              xlim=c(5e2,1e7),
              ylim=c(min(SED$Photom)/3,max(SED$FinalFlux$flux)),
              xlab='Wave / Ang',
              ylab='Flux / Jansky',
              type='l',
              lwd=5,
              col='grey',
              grid=TRUE)
      points(cenwave$cenwave, SED$Photom, pch=16, cex=2, col=rev(rainbow(20,end=2/3)))
      legend('topleft',legend=filters,col=rev(rainbow(20,end=2/3)), pch=16, pt.cex=2)
      if(!is.null(flux)){
        points(x=flux$cenwave, y=flux$flux)
        if(!is.null(flux$fluxerr)){
          magerr(x=flux$cenwave, y=flux$flux, ylo=flux$fluxerr)
        }
      }
   })
   
   output$SED_lum_plot <- renderPlot({
     
     SED=ProSpectSED(SFH=SFHp5,
                  z=0,
                  tau_birth=input$tau_birth,
                  tau_screen=input$tau_screen,
                  tau_AGN=input$tau_AGN,
                  alpha_SF_birth=input$alpha_SF_birth,
                  alpha_SF_screen=input$alpha_SF_screen,
                  alpha_SF_AGN=input$alpha_SF_AGN,
                  AGNlum=10^input$AGNlum,
                  speclib=BC03lr,
                  Dale=Dale_NormTot,
                  AGN=AGN_UnOb_Sparse,
                  burstmass=10^input$burstmass,
                  youngmass=10^input$youngmass,
                  midmass=10^input$midmass,
                  oldmass=10^input$oldmass,
                  ancientmass=10^input$ancientmass
                  )
     
      par(mar=c(3.1,3.1,1.1,1.1))
     
      magicaxis::magplot(SED$FinalLum,
              log='xy',
              xlim=c(5e2,1e7),
              ylim=c(min(SED$FinalLum$lum)*1e3,max(SED$StarsUnAtten$lum, SED$FinalLum$lum)),
              xlab='Wave / Ang',
              ylab='Lum / Lsol/A',
              type='l',
              lwd=5,
              col='grey',
              grid=TRUE)
      lines(SED$StarsUnAtten, col='blue', lwd=2, lty=2)
      lines(SED$StarsAtten, col='green', lwd=2)
      lines(SED$DustEmit, col='darkgreen', lwd=2)
      lines(SED$AGN, col='brown', lwd=2)
      legend('topright', legend=c('Observed Total', 'Unattenuated Stars', 'Attenuated Stars', 'Re-emitted dust', 'AGN'), col=c('grey', 'blue', 'green', 'darkgreen', 'brown'), lty=c(1,2,1,1,1), lwd=c(5,2,2,2,2))
   })
   
   output$SFH_plot <- renderPlot({
     
      burst_SFR=(10^input$burstmass)/1e8
      young_SFR=(10^input$youngmass)/9e8
      mid_SFR=(10^input$midmass)/4e9
      old_SFR=(10^input$oldmass)/4e9
      ancient_SFR=(10^input$ancientmass)/4e9
      
      par(mar=c(3.1,3.1,1.1,1.1))
     
      magicaxis::magplot(c(1, 1, 1e8, 1e8, 1e8, 1e8, 1e9, 1e9, 1e9, 1e9, 5e9, 5e9, 5e9, 5e9,  9e9, 9e9, 9e9, 9e9, 1.3e10, 1.3e10)/1e9,
              c(0, burst_SFR, burst_SFR, 0, 0, young_SFR, young_SFR, 0, 0, mid_SFR, mid_SFR, 0, 0, old_SFR, old_SFR, 0, 0, ancient_SFR, ancient_SFR, 0),
              xlim=c(0,13),
              grid=TRUE,
              xlab='Light Travel Time / Gyr',
              ylab='SFR / Msol/Yr',
              type='s')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
