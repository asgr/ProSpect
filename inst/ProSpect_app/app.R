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
Z=.checkoption('Z',5)

# Define UI for application that draws a histogram
ui = fluidPage(
  # Application title
  fluidRow(
        h1("ProSpect SED")
    ),
  
  # Sidebar with a slider input for number of bins 
  column(2,
         numericInput("z",
                      "Redshift",
                      value = z,
                      min = 0,
                      max = Inf,
                      step=0.1
                      
         ),
         h4('SFH/AGN:'),
         sliderInput("burstmass",
                     "Burst Mass (log10):",
                     min = 1,
                     max = 12,
                     value = burstmass, step=0.1),
         sliderInput("youngmass",
                     "Young Mass (log10):",
                     min = 1,
                     max = 12,
                     value = youngmass, step=0.1),
         sliderInput("midmass",
                     "Mid Mass (log10):",
                     min = 1,
                     max = 12,
                     value = midmass, step=0.1),
         sliderInput("oldmass",
                     "Old Mass (log10):",
                     min = 1,
                     max = 12,
                     value = oldmass, step=0.1),
         sliderInput("ancientmass",
                     "Ancient Mass (log10):",
                     min = 1,
                     max = 12,
                     value = ancientmass, step=0.1),
         sliderInput("AGNlum",
                     "AGN Luminosity (log10):",
                     min = 30,
                     max = 50,
                     value = AGNlum, step=0.2),
         selectInput('Z',
                     'Z',
                     list('1e-04'=1, '4e-04'=2, '0.004'=3, '0.008'=4, '0.02'=5, '0.05'=6),
                     Z
         )
  ),
  column(2,
         h4('Stop and return values:'),
         actionButton("stop","Stop!"),
         h4('Dust inputs:'),
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
                                br(), 
                                plotOutput("SFH_plot", height="250px"),
                                br(),
                                h4('Log10(Mass/Msol):'),
                                DT::dataTableOutput('table')
                              )
                     ),
                     tabPanel('Info',
                              mainPanel(
                                h3("Web Tool Parameter Details"),
                                p(strong("Burst Mass:"),"Total stellar mass formed in a recent burst (0-100 Myrs) in units Msol"),
                                p(strong("Young Mass:"),"Total stellar mass formed in young stars (0.1-1 Gyrs) in units Msol"),
                                p(strong("Mid Mass:"),"Total stellar mass formed in middle aged stars (1-5 Gyrs) in units Msol"),
                                p(strong("Old Mass:"),"Total stellar mass formed in old stars (5-9 Gyrs) in units Msol"),
                                p(strong("Ancient Mass:"),"Total stellar mass formed in ancient stars (9-13 Gyrs) in units Msol"),
                                p(strong("AGN Luminosity:"),"AGN bolometric luminosity in erg/s"),
                                p(strong("Tau Birth"),"Charlot and Fall dust attenuation tau for birth clouds"),
                                p(strong("Tau Screen"),"Charlot and Fall dust attenuation tau for dust screen"),
                                p(strong("Tau AGN"),"Charlot and Fall dust attenuation tau for the AGN torus"),
                                p(strong('Alpha SF Birth'),"Dale dust radiation power law for the birth cloud (lower values mean hotter)"),
                                p(strong('Alpha SF Screen'),"Dale dust radiation power law for the dust screen (lower values mean hotter)"),
                                p(strong('Alpha SF AGN'),"Dale dust radiation power law for the AGN torus (lower values mean hotter)"),
                                p(strong('Z'),"Metallicity of the stars in Fe/H (0.02 is solar)"),
                                
                                h3("Getting Started"),
                                p("New to Synthetic Spectra? Then take a look here: [sedfitting.org](http://www.sedfitting.org/Models.html). There's a bunch of useful stuff there, and a lot of dead links :-(. But it is worth exploring in some detail. In particular for the stellar population models have a look at BASTI, BPASS, Galaxev (which is BC03 to most astronomers), MILES, Pegase, SLUG and Starburst99. For dust stuff have a look at Dale+Helou and Draine+Li, for UV-FIR Da Cunha (which is MagPhys), CIGALE and Grasil. Those are the most popular variants of what they do in their respective fields. Caveats abound about which is better, but these days they are all pretty sophisticated in their own way."),
                                
                                h3("A New ProSpect"),
                                p("ProSpect is a package that aims to help users explore star formation histories (SFH) and spectral energy distributions (SED). Beneath it all it makes use of the BC03 and EMILES sythetic stellar population libraries. On top of this it uses the Charlot and Fall model for birth cloud and screen dust attenuation and re-emits MIR to FIR flux via the incorporation of the most recent Dale dust templates that model the heating of dust by a radiation field and AGN. ProSpect can handle general SFHs via arbitrary functional forms (SFHfunc). The main fixed things to be aware of are the IMF (Chabrier, 2003), SSPs (BC03), dust attention (Charlot and Fall 2000) and dust templates (Dale+ 2014)."),
                                p("Conceptually this probably all sounds a bit similar to MagPhys and Cegale, bar the different dust model (Dale rather than greybody) and some additional SSP libraries. To a degree this is true, but the main reason for putting it all together is to allow proper generative creation of SEDs via the alteration of user accessible parameters and arbitrary SFHs. MagPhys and similar codes do not allow easy access to the under-the-hood generative functionality which might allow this. Doing this is with a longer term aim of incorporating ProSpect into ProFit for multiband morphological decomposition of galaxies, revealing their component-wise SFHs. Before the completion of that component of the project ProSpect offers a readily accessible interface to multiband fitting of SED in order to measure (e.g.) stellar and dust masses, and it can also applied to arbitary SFHs computed by other codes, e.g. create a realistic SED for a semi-analytic (SAM) model SFH."),
                                p("
The 5 phase model (SFHp5) covers 5 key phases of star formation and models them as having constant star formation in each of these periods. By default they cover 0 - 100 Myr (burst, generally considered the shortest phase that broad band photometry can be sensitive to); 0.1-1 Gyr (young; dominated by hot young stars and violent phases of stellar evolution); 1-5 Gyr (mid), 5-9 Gyr (old) and 9-13 Gyr (ancient). It is not really possible to break the SFH into any more independent phases than this using broad band photometry alone, but a physically motivated functional model (e.g. an exponentially decling SFR, or constant SFR) can be used using the functional interface in SFHfunc.")
                              )
                     )
         )
  )
)

server = function(input, output) {
  
  observe({
    if(input$stop > 0){
      stopApp(returnValue = list(
        z=max(input$z,2.261565e-09,na.rm=TRUE),
        burstmass=10^input$burstmass,
        youngmass=10^input$youngmass,
        midmass=10^input$midmass,
        oldmass=10^input$oldmass,
        ancientmass=10^input$ancientmass,
        AGNlum=10^input$AGNlum,
        tau_birth=input$tau_birth,
        tau_screen=input$tau_screen,
        tau_AGN=input$tau_AGN,
        alpha_SF_birth=input$alpha_SF_birth,
        alpha_SF_screen=input$alpha_SF_screen,
        alpha_SF_AGN=input$alpha_SF_AGN,
        Z=as.integer(input$Z)
      )
      )
    }
  })
  
  output$SED_flux_plot = renderPlot({
    
    SED=ProSpectSED(SFH=SFHfunc,
                    z=max(input$z,2.261565e-09,na.rm=TRUE),
                    m1=10^input$burstmass/1e8,
                    m2=10^input$youngmass/9e8,
                    m3=10^input$midmass/4e9,
                    m4=10^input$oldmass/4e9,
                    m5=10^input$ancientmass/4e9,
                    AGNlum=10^input$AGNlum,
                    tau_birth=input$tau_birth,
                    tau_screen=input$tau_screen,
                    tau_AGN=input$tau_AGN,
                    alpha_SF_birth=input$alpha_SF_birth,
                    alpha_SF_screen=input$alpha_SF_screen,
                    alpha_SF_AGN=input$alpha_SF_AGN,
                    Z=as.integer(input$Z),
                    speclib=BC03lr,
                    Dale=Dale_NormTot,
                    AGN=AGN_UnOb_Sparse,
                    filtout=filtout
    )
    
    par(mar=c(3.1,3.1,1.1,1.1))
    
    magicaxis::magplot(SED$FinalFlux,
                       log='xy',
                       xlim=c(5e2,1e7),
                       ylim=c(min(SED$Photom)/3,max(SED$FinalFlux$flux)),
                       xlab='Wavelength (Angstrom)',
                       ylab='Flux Density (Jansky)',
                       type='l',
                       lwd=5,
                       col='grey',
                       grid=TRUE)
    colvec=rev(rainbow(20,s=c(0.5,1),v=seq(0.5,1,len=4),end=5/6))
    points(cenwave$cenwave, SED$Photom, pch=16, cex=2, col=colvec)
    legend('topleft',legend=filters,col=colvec, pch=16, pt.cex=2)
    if(!is.null(flux)){
      points(x=flux$cenwave, y=flux$flux)
      if(!is.null(flux$fluxerr)){
        magerr(x=flux$cenwave, y=flux$flux, ylo=flux$fluxerr)
      }
    }
  })
  
  output$SED_lum_plot = renderPlot({
    
    SED=ProSpectSED(SFH=SFHfunc,
                    z=max(input$z,2.261565e-09,na.rm=TRUE),
                    m1=10^input$burstmass/1e8,
                    m2=10^input$youngmass/9e8,
                    m3=10^input$midmass/4e9,
                    m4=10^input$oldmass/4e9,
                    m5=10^input$ancientmass/4e9,
                    AGNlum=10^input$AGNlum,
                    tau_birth=input$tau_birth,
                    tau_screen=input$tau_screen,
                    tau_AGN=input$tau_AGN,
                    alpha_SF_birth=input$alpha_SF_birth,
                    alpha_SF_screen=input$alpha_SF_screen,
                    alpha_SF_AGN=input$alpha_SF_AGN,
                    Z=as.integer(input$Z),
                    speclib=BC03lr,
                    Dale=Dale_NormTot,
                    AGN=AGN_UnOb_Sparse,
                    filtout=filtout
    )
    
    par(mar=c(3.1,3.1,1.1,1.1))
    
    magicaxis::magplot(SED$FinalLum,
                       log='xy',
                       xlim=c(5e2,1e7),
                       ylim=c(min(SED$FinalLum$lum)*1e3,max(SED$StarsUnAtten$lum, SED$FinalLum$lum)),
                      xlab='Wavelength (Angstrom)',
                       ylab='Luminosity Density (Lsol/A)',
                       type='l',
                       lwd=5,
                       col='grey',
                       grid=TRUE)
    lines(SED$StarsUnAtten, col='blue', lwd=2, lty=2)
    lines(SED$StarsAtten, col='green', lwd=2)
    lines(SED$DustEmit, col='darkgreen', lwd=2)
    lines(SED$AGN, col='brown', lwd=2)
    legend('topright', legend=c('Observed Total', 'Unattenuated Stars', 'Attenuated Stars', 'Re-emitted dust', 'AGN'), col=c('grey', 'blue', 'green', 'darkgreen', 'brown'), lty=c(1,2,1,1,1), lwd=c(5,2,2,2,2))
    
    SMphases=SMstarfunc(z=max(input$z,2.261565e-09,na.rm=TRUE),
                        m1=10^input$burstmass/1e8,
                        m2=10^input$youngmass/9e8,
                        m3=10^input$midmass/4e9,
                        m4=10^input$oldmass/4e9,
                        m5=10^input$ancientmass/4e9,
                        Z=as.integer(input$Z)
    )
    
    SMphases=round(log10(SMphases),2)
    
    SMphases=data.frame(c(SMphases['TotSMform'],SMphases[1:5]),c(SMphases['TotSMstar'],SMphases[1:5+5]))
    SMphases=cbind(c('Total','Burst','Young','Mid','Old','Ancient'),SMphases)
    colnames(SMphases)=c('Phase','Formed','Remaining')
    
    output$table = DT::renderDataTable(SMphases)
  })
  
  output$SFH_plot = renderPlot({
    
    TravelTime=cosdistTravelTime(z=max(input$z,2.261565e-09,na.rm=TRUE), H0 = 67.8, OmegaM = 0.308)
    
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
                       xlab='Light Travel Time (Gyr)',
                       ylab='SFR (Msol/Yr)',
                       type='s')
    abline(v=13.8-TravelTime, lty=2)
    rect(13.8-TravelTime, 0, 100, 100, col=hsv(alpha=0.2))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
