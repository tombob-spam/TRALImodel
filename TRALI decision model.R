#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(tidyverse)
library(reshape2)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("TRALI decision model"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fluidRow("Parameters"),
            fluidRow(sliderInput("pd",
                        "Donor cognate Ab prevalence",
                        min = 0,
                        max = 1,
                        value = 0.05),
            sliderInput("pp",
                        "P(Ab+|Antibody mediated TRALI)",
                        min = 0,
                        max = 1,
                        value = 0.9),
            numericInput("pt", "P(recurrent TRALI)", 0.001, min=0, max=1)),
            
          
            fluidRow("Utility values"),
            
            fluidRow(
              column(6, numericInput("u0", "Temporary deferral of 1 donor", 1, min=1)
                     
              ),
              column(6, numericInput("ud", "Permanent deferral of 1 donor", 5, min=0)
                     
              )
            ),
            fluidRow(
              column(6, numericInput("uc", "Diagnostic", 5, min=0)
                     
              ),
              column(6, numericInput("ut", "Preventable TRALI", 10000, min=0)
                     
              )
            )
           
        ),
       
       
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           verbatimTextOutput("capture")
        )
    )
)
#Define functions for decision model
UtilityGain<-function(p, n, pd, pp, pt, ud, ut, uc){
  uc * pp * (1 - p) * p * (1 - pd) ^ n + ut * (1 - (1 - pt) ^ (p + (n - p) * pd)) - ud * (p - pd * (p - n)) - n 
}

UtilityDeferAll<-function(p, n, pd, pp, pt, ud, ut, uc){
  (1 - (1 - pt) ^ (p + pd * (n - p))) * ut - n * ud
}

Decision<-function(p, n, pd, pp, pt, ud, ut, uc){
  sign_T<- sign(UtilityGain(p, n, pd, pp, pt, ud, ut, uc) - UtilityDeferAll(p, n, pd, pp, pt, ud, ut, uc));
  pos_T = max(UtilityGain(p, n, pd, pp, pt, ud, ut, uc), UtilityDeferAll(p, n, pd, pp, pt, ud, ut, uc)) > 0;
  return(sign_T * pos_T)
}

DefaultDecision<-function(p,n){Decision(p,n,0.05,0.9,0.001,5,10000,5)}




# Define server logic required to draw a histogram
server <- function(input, output) {
  
    DecisionTable<-reactive({
      CurryDecision<-function(p,n){Decision(p,n,input$pd,input$pp,input$pt,input$ud/input$u0,input$ut/input$u0,input$uc/input$u0)}
      
      p<-seq(0,1,by=0.1)
      n<-1:20
      z<-outer(p,n,FUN=Vectorize(CurryDecision))
      z<-z[11:1,]
      colnames(z)<-n
      rownames(z)<-p[11:1]
      
      zz<-melt(z)
      zz$value<-as.factor(zz$value)
      zz
      })
    
    
  
    output$distPlot <- renderPlot(
      ggplot(data = DecisionTable(), aes(x=Var2, y=Var1, fill=value)) + 
        geom_tile()+
        ggtitle("Recommended action")+
        xlab("number of donors")+
        ylab("pretest probability")+
        
        scale_fill_manual(values=c("black", "yellow", "red"), 
                          name="Action",
                          breaks=c(-1, 0, 1),
                          labels=c("Defer", "No test", "Test"))
    )
    output$capture=renderText("This model must not be used for clinical decision making \nwithout agreeing locally appropriate parameters and utility values ")
  

}
# Run the application 
shinyApp(ui = ui, server = server)
