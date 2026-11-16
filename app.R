#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(readr)
library(flan)
library(ggplot2)
library(ggforce)

#read_csv("AppDetails2.csv")->table
#table[,1:4]->table
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("BIOL336: Mutation Rate Starter"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            shinyWidgets::sliderTextInput("m","Mutation Rate:",
                                          choices=c( 1.00e-07, 2e-07, 
                                                     5e-07 ,8e-07, 1e-06 ,
                                                     2e-06, 4e-06 ,6e-06,
                                                    1.00e-05),
                                          selected= 1e-06, grid = T),
            sliderInput(inputId = "s", label = "Mutant Fitness", value = 1, 
                      min = 0.2, max = 1.5, step = 0.1)#,
           
            # shinyWidgets::sliderTextInput("d","Death Rate:",
            #                               choices=c( 0,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1),
            #                               selected= 0, grid = T)#  ,
            # actionButton("update" ,"Update View", icon("refresh"),
            #              class = "btn btn-primary")
            ,
            tags$figure(
              class = "centerFigure",
              tags$img(
                src = "diagram.png",
                width = 400,
                alt = "Cartoon of population growth during which mutations accumulate followed by selection for mutants"
              ),
              tags$figcaption("")
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
          #plotOutput("plot_phylo"),
          plotOutput("plot"),
           #plotOutput("distPlot")
          #textOutput("text"),
          verbatimTextOutput("verb"),
          verbatimTextOutput("verb2")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  
  # output$text <- renderText({ 
  #   rflan(1e4,mutprob=input$m,mfn=1e6,
  #         fitness = 1/input$s,death=input$d,
  #         dist=list("lnorm",meanlog=-0.3795851,sdlog=0.3016223))$mc->mc
  #   paste0("Median Number of Survivors: ",median(mc)) })
  output$verb <- renderText({ 
    rflan(1e4,mutprob=input$m,mfn=1e6,
          fitness = 1/input$s,#death=input$d,
          dist=list("lnorm",meanlog=-0.3795851,sdlog=0.3016223))$mc->mc
    paste0("Median Number of Survivors: ",median(mc)," (N trials = 10,000)"  )})
  output$verb2 <- renderText({ 
    dflan(0,mutations=input$m*1e6,
          fitness = 1/input$s)%>%signif(3)->p0
    paste0("Probability of ≥1 survivor: ",1-p0)  })
  
  # output$plot_phylo<-renderPlot({
  #   flan::draw.clone(t=6,col=c('black','red'),mutprob=0.06,
  #                    dist=list(name="dirac",location=1))
  # })
  
  output$plot<-renderPlot({
    rflan(18,mutprob=input$m,mfn=1e6,
                     fitness = 1/input$s,#death=input$d,
                     dist=list("lnorm",meanlog=-0.3795851,sdlog=0.3016223))$mc->mc
    
    
    ## set up function
    vis_points_inside_circle <- function(n) {
      
      # part 1 -- set up empty circle
      df_empty_circle <-
        data.frame(x = 0,
                   y = 0, 
                   r = 1)
      
      p_empty_circle <-
        ggplot(df_empty_circle) +
        geom_circle(aes(x0 = x, y0 = y, r = r)) + 
        coord_fixed() +
        theme_void()
      
      # part 2 -- set up points scatter
      r <- runif(n)
      th <- runif(n, 0, 2 * pi)
      
      df_circular_points_scatter <- 
        data.frame(x = sqrt(r) * cos(th), 
                   y = sqrt(r) * sin(th)
                   
        )
      
      
      # part 3 -- combine circle and points
      p_empty_circle +
        geom_point(data = df_circular_points_scatter, 
                   aes(x = x, y = y),
                   col='red')+
        ggtitle(n)
    }
    
    ## visualize
    library(gridExtra)
    
    set.seed(2021)
    p1000_1  <- vis_points_inside_circle(n = mc[1 ])
    p1000_2  <- vis_points_inside_circle(n = mc[2 ])
    p1000_3  <- vis_points_inside_circle(n = mc[3 ])
    p1000_4  <- vis_points_inside_circle(n = mc[4 ])
    p1000_5  <- vis_points_inside_circle(n = mc[5 ])
    p1000_6  <- vis_points_inside_circle(n = mc[6 ])
    p1000_7  <- vis_points_inside_circle(n = mc[7 ])
    p1000_8  <- vis_points_inside_circle(n = mc[8 ])
    p1000_9  <- vis_points_inside_circle(n = mc[9 ])
    p1000_10 <- vis_points_inside_circle(n = mc[10])
    p1000_11 <- vis_points_inside_circle(n = mc[11])
    p1000_12 <- vis_points_inside_circle(n = mc[12])
    p1000_13 <- vis_points_inside_circle(n = mc[13])
    p1000_14 <- vis_points_inside_circle(n = mc[14])
    p1000_15 <- vis_points_inside_circle(n = mc[15])
    p1000_16 <- vis_points_inside_circle(n = mc[16])
    p1000_17 <- vis_points_inside_circle(n = mc[17])
    p1000_18 <- vis_points_inside_circle(n = mc[18])
    
    
    gridExtra::grid.arrange(p1000_1 ,
                            p1000_2 ,
                            p1000_3 ,
                            p1000_4 ,
                            p1000_5 ,
                            p1000_6 ,
                            p1000_7 ,
                            p1000_8 ,
                            p1000_9 ,
                            p1000_10,
                            p1000_11,
                            p1000_12, 
                            p1000_13,
                            p1000_14,
                            p1000_15,
                            p1000_16,
                            p1000_17,
                            p1000_18,
                            nrow = 3,
                            top="18 Example Populations (1 million individuals)\n Mutant Survivors Shown")
    
    
  })

    # output$distPlot <- renderPlot({
    #     # generate bins based on input$bins from ui.R
    #     x    <- faithful[, 2]
    #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
    # 
    #     # draw the histogram with the specified number of bins
    #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
    #          xlab = 'Waiting time to next eruption (in mins)',
    #          main = 'Histogram of waiting times')
    # })
  
  # Source - https://stackoverflow.com/q/68619095
  # Posted by Emman, modified by community. See post 'Timeline' for change history
  # Retrieved 2026-11-14, License - CC BY-SA 4.0
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

