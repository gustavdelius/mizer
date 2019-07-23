library(tidyverse)
library(ggplot2)
library(shiny)
library(miniUI)

params <- readRDS("inst/humboldt/params0707.rds")
alpha <- 1/3

ppmr <- read_csv("inst/humboldt/Data_PPMR.csv") %>% 
    mutate(Species = str_to_title(Species),
           logpredprey = log(wpredator / wprey)) %>% 
    group_by(Species) %>% 
    mutate(weight = Nprey / sum(Nprey),
           weight_biomass = Nprey * wprey / sum(Nprey * wprey),
           weight_kernel = Nprey / wprey^(1 + alpha - params@lambda),
           weight_kernel = weight_kernel / sum(weight_kernel))

sigmoid <- function(x, p, s) {
    1/(1 + exp((p - x) / s))
}

double_sigmoid <- function(x, p_l, s_l, p_r, s_r, ex) {
    sigmoid(x, p_l, s_l) *
        exp(-ex * x) *
        sigmoid(-x, -p_r, s_r)
}

dx <- log(params@w[2]/params@w[1])
adjust <- 1/8

fitKernel <- function(ppmr, species) {
    ppmr <- ppmr %>% 
        filter(Species == species)
    
    p_l <- min(ppmr$logpredprey)
    p_r <- max(ppmr$logpredprey)
    s_l <- 0.5
    s_r <- 0.5
    ex <- 0.1
    
    x <- seq(0, max(ppmr$logpredprey) * 1.1, by = dx)
    
    ui <- fluidPage(sidebarLayout(
        
        ## Sidebar ####
        sidebarPanel(
            actionButton("done", "Done",
                         onclick = "setTimeout(function(){window.close();},500);"),
            sliderInput("p_l", "p_l",
                        value = p_l,
                        min = p_l - 1,
                        max = p_l + 5),
            sliderInput("s_l", "s_l",
                        value = s_l,
                        min = 0.01,
                        max = 2),
            sliderInput("p_r", "p_r",
                        value = p_r,
                        min = p_r - 5,
                        max = p_r + 1),
            sliderInput("s_r", "s_r",
                        value = s_r,
                        min = 0.01,
                        max = 2),
            sliderInput("ex", "ex",
                        value = ex,
                        min = 0,
                        max = 1)
        ),
        
        mainPanel(
            plotOutput("plot"),
            sliderInput("adjust", "Adjust bandwidth",
                        value = 1,
                        min = 0,
                        max = 2,
                        step = 0.1)
        )
    ))
    
    server <- function(input, output, session) {
        
        # Render the plot
        output$plot <- renderPlot({
            pl <- ggplot(ppmr) +
                geom_density(aes(logpredprey, weight = weight,
                                 colour = "Numbers"),
                             adjust = input$adjust) +
                geom_density(aes(logpredprey, weight = weight_biomass,
                                 colour = "Biomass"),
                             adjust = input$adjust) +
                # geom_density(aes(logpredprey, weight = weight_kernel,
                #                  colour = "Kernel"),
                #              adjust = adjust) +
                xlab("Log of predator/prey mass ratio") +
                ylab("Density") +
                scale_colour_manual(values = c(Numbers = "Blue",
                                               Biomass = "Green",
                                               Kernel = "Black")) +
                expand_limits(x = 0)
            df <- tibble(
                x = x,
                Kernel = double_sigmoid(
                    x, 
                    p_l = input$p_l, 
                    s_l = input$s_l, 
                    p_r = input$p_r, 
                    s_r = input$s_r, 
                    ex = input$ex)) %>% 
                mutate(Numbers = Kernel / exp((1 + alpha - params@lambda) * x),
                       Biomass = Numbers / exp(x),
                       Kernel = Kernel / sum(Kernel) / dx,
                       Numbers = Numbers / sum(Numbers) / dx,
                       Biomass = Biomass / sum(Biomass) / dx) %>% 
                gather(key = Type, value = "Density",
                       Numbers, Biomass)
            
            pl + geom_line(data = df,
                           aes(x, Density, colour = Type),
                           size = 3) 
        })
        
        # Handle the Done button being pressed.
        observeEvent(input$done, {
            stopApp(list(
                p_l = input$p_l, 
                s_l = input$s_l, 
                p_r = input$p_r, 
                s_r = input$s_r, 
                ex = input$ex))
        })
    }
    
    runGadget(ui, server, viewer = browserViewer())
}