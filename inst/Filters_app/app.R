#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ProSpect)
library(magicaxis)
data("ProFiltTrans_Shark")

filter_names = ProFiltTrans_Shark$Names

z_steps = round((ProFiltTrans_Shark$zsteps[2:length(ProFiltTrans_Shark$zsteps)] + ProFiltTrans_Shark$zsteps[1:(length(ProFiltTrans_Shark$zsteps)-1)])/2,2)
z_choices = lapply(1:length(z_steps),function(x){x})
names(z_choices) = z_steps

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("ProSpect Filter Transforms"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("filter_target",
                        "Target Filter:",
                        choices = filter_names,
                        selected = 'r_SDSS'
                        ),
            selectInput("filter_ref1",
                        "Reference Filter 1:",
                        choices = filter_names,
                        selected = 'g_VST'
            ),
            selectInput("filter_ref2",
                        "Reference Filter 2:",
                        choices = filter_names,
                        selected = 'r_VST'
            ),
            selectInput("z_step",
                        "Redshift:",
                        choices = z_choices,
                        selected = 1
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            uiOutput("trans_func"),
            plotOutput("filter_comp"),
            plotOutput("trans_hist")
            #plotOutput("trans_plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    transform=reactive({
        filterTranMags(
            mag_in = ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,c(input$filter_ref1, input$filter_ref2)],
            mag_out = ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target],
            return = 'bestall'
        )
    })
    
    ref1_offset=reactive({
        temp=ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target] - ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_ref1]
        return(c(med=median(temp,na.rm = TRUE), sd=sd(temp)))
    })
    
    ref2_offset=reactive({
        temp=ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target] - ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_ref2]
        return(c(med=median(temp,na.rm = TRUE), sd=sd(temp)))
    })
    
    output$trans_func = renderUI({
        list(
            HTML(input$filter_target,'~',names(transform()$params[1]),'<br/>'),
            HTML('alpha =', round(as.numeric(transform()$params[[1]]['alpha']), digits=4), '<br/>'),
            HTML('beta =', round(as.numeric(transform()$params[[1]]['beta']), digits=4), '<br/>'),
            HTML('sigma =', round(as.numeric(transform()$params[[1]]['sigma']), digits=4), '<br/><br/>'),
            HTML(input$filter_target,'~', input$filter_ref1, switch(sign(ref1_offset()[1])+2,'-',NA,'+'), round(abs(as.numeric(ref1_offset()[1])), digits=4), '+/-', round(as.numeric(ref1_offset()[2]), digits=4), '<br/>'),
            HTML(input$filter_target,'~', input$filter_ref2, switch(sign(ref2_offset()[1])+2,'-',NA,'+'), round(abs(as.numeric(ref2_offset()[1])), digits=4), '+/-', round(as.numeric(ref2_offset()[2]), digits=4), '<br/>')
        )
    })
    
    output$filter_comp = renderPlot({
        if(input$filter_target %in% pivwave$filter){
            filttemp = getfilt(input$filter_target)    
        } else{
            filttemp = EAZY_filters$filters[[which(EAZY_filters$info == input$filter_target)]]
        }
        filttemp[,2] = filttemp[,2] / max(filttemp[,2], na.rm=TRUE)
        magplot(
            filttemp,
            type = 'l',
            xlim = c(2e3,3e4),
            ylim = c(0,1),
            log = 'x',
            grid = TRUE,
            xlab = 'Wavelength / Ang',
            ylab = 'Response'
        )
        
        if(input$filter_ref1 %in% pivwave$filter){
            filttemp = getfilt(input$filter_ref1)    
        } else{
            filttemp = EAZY_filters$filters[[which(EAZY_filters$info == input$filter_ref1)]]
        }
        filttemp[,2] = filttemp[,2] / max(filttemp[,2], na.rm=TRUE)
        lines(
            filttemp,
            col = 'blue'
        )
        
        if(input$filter_ref2 %in% pivwave$filter){
            filttemp = getfilt(input$filter_ref2)    
        } else{
            filttemp = EAZY_filters$filters[[which(EAZY_filters$info == input$filter_ref2)]]
        }
        filttemp[,2] = filttemp[,2] / max(filttemp[,2], na.rm=TRUE)
        lines(
            filttemp,
            col = 'red'
        )
        legend(
            'topleft',
            legend = c(input$filter_target, input$filter_ref1, input$filter_ref2),
            col = c('black', 'blue', 'red'),
            lty = 1
        )
    })
    
    output$trans_hist = renderPlot({
        magplot(
            density(ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target] - transform()$predict),
            xlab = paste0(input$filter_target,' - Transformed / Ref1 / Ref2'),
            ylab = 'PDF',
            xlim = c(-0.1,0.1),
            grid = TRUE
        )
        lines(
            density(ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target] - ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_ref1]),
            col = 'blue'
        )
        abline(v=ref1_offset()[1], col='blue')
        abline(v=ref1_offset()[1] - ref1_offset()[2], col='blue', lty=3)
        abline(v=ref1_offset()[1] + ref1_offset()[2], col='blue', lty=3)
        lines(
            density(ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target] - ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_ref2]),
            col = 'red'
        )
        abline(v=ref2_offset()[1], col='red')
        abline(v=ref2_offset()[1] - ref2_offset()[2], col='red', lty=3)
        abline(v=ref2_offset()[1] + ref2_offset()[2], col='red', lty=3)
        legend(
            'topleft',
            legend = c('Transformed', input$filter_ref1, input$filter_ref2),
            col = c('black', 'blue', 'red'),
            lty = 1
        )
    })
    
    # output$trans_plot = renderPlot({
    #     magplot(
    #         ProFiltTrans_Shark$maglist[[as.integer(input$z_step)]][,input$filter_target],
    #         transform()$predict,
    #         xlab=input$filter_target,
    #         ylab='Prediction',
    #         grid=TRUE
    #     )
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
