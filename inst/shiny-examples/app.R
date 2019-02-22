BC03lr=NULL
data('BC03lr', envir = environment())
Dale_NormTot=NULL
data('Dale_NormTot', envir = environment())
AGN_UnOb_Sparse=NULL
data('AGN_UnOb_Sparse', envir = environment())

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("ProSpect SED"),
   
   # Sidebar with a slider input for number of bins 
        column(2,
         sliderInput("burstmass",
                     "Burst Mass:",
                     min = 1,
                     max = 12,
                     value = 8, step=0.1),
         sliderInput("youngmass",
                     "Young Mass:",
                     min = 1,
                     max = 12,
                     value = 9, step=0.1),
         sliderInput("midmass",
                     "Mid Mass:",
                     min = 1,
                     max = 12,
                     value = 9, step=0.1),
         sliderInput("oldmass",
                     "Old Mass:",
                     min = 1,
                     max = 12,
                     value = 9, step=0.1),
         sliderInput("ancientmass",
                     "Ancient Mass:",
                     min = 1,
                     max = 12,
                     value = 9, step=0.1),
         sliderInput("AGNlum",
                     "AGN Luminosity:",
                     min = 40,
                     max = 47,
                     value = 42, step=0.1)
        ),
        column(2,
         sliderInput("tau_birth",
                     "Tau Birth:",
                     min = 0,
                     max = 2,
                     value = 1, step=0.1),
         sliderInput("tau_screen",
                     "Tau Screen:",
                     min = 0,
                     max = 2,
                     value = 0.3, step=0.1),
         sliderInput("tau_AGN",
                     "Tau AGN:",
                     min = 0,
                     max = 2,
                     value = 1, step=0.1),
         sliderInput("alpha_SF_birth",
                     "Alpha SF Birth:",
                     min = 0,
                     max = 4,
                     value = 1, step=0.1),
         sliderInput("alpha_SF_screen",
                     "Alpha SF Screen:",
                     min = 0,
                     max = 4,
                     value = 3, step=0.1),
         sliderInput("alpha_SF_AGN",
                     "Alpha SF AGN:",
                     min = 0,
                     max = 4,
                     value = 0, step=0.1)
        ),
      
      # Show a plot of the generated distribution
      column(8, 
        mainPanel(
           plotOutput("SED_plot", height="600px")
        )
      )
)

server <- function(input, output) {
  
   output$SED_plot <- renderPlot({
     
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
     
      magicaxis::magplot(SED$FinalLum,
              log='xy',
              xlim=c(5e2,1e7),
              ylim=c(5e2,1e7),
              xlab='Wave / Ang',
              ylab='Lum / Lsol/A',
              type='l',
              lwd=5)
      lines(SED$StarsUnAtten, col='blue', lty=2)
      lines(SED$StarsAtten, col='green')
      lines(SED$DustEmit, col='darkgreen')
      lines(SED$AGN, col='brown')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
