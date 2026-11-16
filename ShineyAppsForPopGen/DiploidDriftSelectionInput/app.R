#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

sim_sel <- function( n_gen = 10,
                     n_reps = 5,
                     p0 = 0.5,
                     N = 100,
                     wAA = 0.6,
                     wAa = 1,
                     waa = 0.2
){
  p_ <- matrix(NA, nrow = n_gen, ncol = n_reps)
  colnames(p_) <- paste0("sim ", 1:n_reps)
  p_[1,] <- rep(p0, n_reps)
  
  all_t <- 2:n_gen
  for(t in all_t){
    for(j in 1:n_reps){
      ave_W <- wAA * (p_[t-1,j])^2 +
        wAa * 2 * p_[t-1,j] * (1 - p_[t-1,j]) +
        waa * (1 - p_[t-1,j])^2
      det_p <- (wAA * (p_[t-1,j])^2 + wAa * p_[t-1,j] * (1 - p_[t-1,j])) / ave_W
      p_[t,j] <- rbinom(n = 1, 2*N, prob = det_p) / (2*N). 
# Note that we have to use 2*N when sampling within a diploid population
    }
  }
  data.frame(
    cbind(
      time = c(1,all_t), 
      pA = p_
    )
  ) %>%
    gather(sim, pA, -time)
}

#user interface
ui <- pageWithSidebar( 
  
  headerPanel = headerPanel("BIOL336: Diploid model of selection with drift"),
  
  sidebarPanel( 
    
    HTML("<p style='font-size:14px'><B>Frequency of allele A over time (p).</B>"),
    HTML("<p style='font-size:10px'>Change settings and then hit Go.</p>"),
    
    actionButton("goButton", "GO"),
    
    numericInput(inputId = "N", label = "Number of individuals", value = 50,
                min = 5, max = 1000, step = 5),
    numericInput(inputId = "n_reps", label = "Number of replicates", value = 10,
                min = 1, max = 100, step = 1),
    numericInput(inputId = "n_gen", label = "Number of generations", value = 50,
                min = 0, max = 100, step = 10),
    numericInput(inputId = "wAA", label = "Fitness of AA", value = 1,
                min = 0, max = 2, step = 0.01),
    numericInput(inputId = "wAa", label = "Fitness of Aa", value = 1,
                min = 0, max = 2, step = 0.01),
    numericInput(inputId = "waa", label = "Fitness of aa", value = 1,
                min = 0, max = 2, step = 0.01),
    numericInput(inputId = "p0", label = "Initial allele frequencies", value = 0.5,
                min = 0, max = 1, step = 0.01),

    downloadButton(outputId = 'sim_out.csv', label = 'Download'),
    
    HTML("<p> </p>"),
    HTML("<p style='font-size:12px'>Alternatively, <A href='https://shiney.zoology.ubc.ca/otto/DiploidDriftSelection'>click here</a> if you would like to input numbers with sliders.<P>"), 
    
    HTML("<p> </p>"),
    HTML("<p style='font-size:12px'>Take me to:
         <UL><LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/HaploidSelection/'>Haploid selection</a> (<A href='https://shiney.zoology.ubc.ca/otto/HaploidMutationSelection/'>with mutation</a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidSelection/'>Diploid selection</a> (<A href='https://shiney.zoology.ubc.ca/otto/DiploidMutationSelection/'>with mutation</a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidDriftSelection/'>Diploid selection with drift</a>
         </UL><P>"),
    HTML("<p style='font-size:8px'>Download R <A href='https://www.zoology.ubc.ca/~otto/Research/ShineyAppsForPopGen.zip'>source code</a>.<BR>
          Modified from: Copyright (c) 2017 Silas Tittes, MIT License, <A href='https://github.com/silastittes/shiny_popgen'>https://github.com/silastittes/shiny_popgen</a></p>")
    
  ), 
  mainPanel =  mainPanel(
    plotOutput(outputId = 'fitness')
  )
)

#back end code and response to user input
server <- function(input, output){
  
  #each time user hits "go"
  rand <- eventReactive(input$goButton, {
    
    sim_sel(input$n_gen, input$n_reps, input$p0, input$N, input$wAA, input$wAa, input$waa)
    
  })
  
  output$fitness <- renderPlot({
    
    out <- rand()
    
    out %>%
      ggplot(aes(time-1, pA, colour = sim)) +
      geom_line() +
      xlab("Generation") +
      ylab("Allele frequency, p") +
      ylim(0, 1) +
      geom_hline(yintercept = out$pA[1], lty = 2) +
#      theme(panel.background = element_rect(
#        fill =  rgb(30, 144, 255, 25, maxColorValue = 255)),
#        text = element_text(size=16, family= "Times"
#        )
#  )
      theme_classic()+theme(text = element_text(size=20))

  })
  
  output$sim_out.csv <- downloadHandler(
    filename = function() {'sim-out.csv'},
    content = function(file) {
      write.csv(rand(), file, row.names = F, quote = F)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)