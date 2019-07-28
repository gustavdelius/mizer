library(tidyverse)
library(ggplot2)
library(shiny)
library(assertthat)

kernel_fn <- power_law_pred_kernel

fitKernel <- function(stomach, species_params = data.frame(), 
                      dx, alpha = 1/3, lambda = 2) {
    
    assert_that(is.data.frame(stomach),
                all(c("wprey", "wpredator") %in% names(stomach)),
                is.data.frame(species_params))
    
    if (!("Nprey" %in% names(stomach))) stomach$Nprey <- 1
    
    stomach <- stomach %>% 
        mutate(logpredprey = log10(wpredator / wprey),
               weight = Nprey / sum(Nprey),
               weight_biomass = Nprey * wprey / sum(Nprey * wprey),
               weight_kernel = Nprey / wprey^(1 + alpha - lambda),
               weight_kernel = weight_kernel / sum(weight_kernel))
    
    vars <- names(species_params)
    species_params <-
        set_species_param_default(species_params, "kernel_l_l", min(stomach$logpredprey))
    species_params <-
        set_species_param_default(species_params, "kernel_l_r", max(stomach$logpredprey))
    species_params <-
        set_species_param_default(species_params, "kernel_u_l", 2)
    species_params <-
        set_species_param_default(species_params, "kernel_u_r", 2)
    species_params <-
        set_species_param_default(species_params, "kernel_exp", -0.5)
    
    x <- seq(0, max(stomach$logpredprey) * 1.1, by = dx)
    r <- 10^x
    
    ui <- fluidPage(sidebarLayout(
        
        ## Sidebar ####
        sidebarPanel(
            actionButton("done", "Done",
                         onclick = "setTimeout(function(){window.close();},500);"),
            sliderInput("l_l", "l_l",
                        value = species_params$kernel_l_l,
                        min = round(species_params$kernel_l_l) - 2,
                        max = round(species_params$kernel_l_l) + 2,
                        step = 0.1),
            sliderInput("u_l", "u_l",
                        value = species_params$kernel_u_l,
                        min = 1,
                        max = 10,
                        step = 0.1),
            sliderInput("l_r", "l_r",
                        value = species_params$kernel_l_r,
                        min = round(species_params$kernel_l_r) - 2,
                        max = round(species_params$kernel_l_r) + 2,
                        step = 0.1),
            sliderInput("u_r", "u_r",
                        value = species_params$kernel_u_r,
                        min = 1,
                        max = 5,
                        step = 0.1),
            sliderInput("ex", "exp",
                        value = species_params$kernel_exp,
                        min = -1,
                        max = 0),
            hr(),
            sliderInput("s", "scaling",
                        value = 1,
                        min = 0,
                        max = 1)
        ),
        
        mainPanel(
            plotOutput("plot"),
            sliderInput("adjust", "Adjust bandwidth",
                        value = 0.5,
                        min = 0.1,
                        max = 1,
                        step = 0.1),
            plotOutput("plot_kernel")
        )
    ))
    
    server <- function(input, output, session) {
        
        # Render the plot
        output$plot <- renderPlot({
            
            stomach <- stomach %>% 
                mutate(weight_scaled = Nprey * wprey^input$s,
                       weight_scaled = weight_scaled / sum(weight_scaled))
            pl <- ggplot(stomach) +
                geom_density(aes(logpredprey, weight = weight_scaled,
                                 colour = "Scaled"),
                             adjust = input$adjust) +
                geom_density(aes(logpredprey, weight = weight,
                                 colour = "Numbers"),
                             adjust = input$adjust) +
                geom_density(aes(logpredprey, weight = weight_biomass,
                                 colour = "Biomass"),
                             adjust = input$adjust) +
                xlab("Log_10 of predator/prey mass ratio") +
                ylab("Density") +
                scale_colour_manual(values = c(Numbers = "Blue",
                                               Biomass = "Green",
                                               Scaled = "Black")) +
                expand_limits(x = 0)
            df <- tibble(
                x = x,
                Kernel = power_law_pred_kernel(
                    r, 
                    kernel_exp = input$ex,
                    kernel_l_l = input$l_l,
                    kernel_u_l = input$u_l,
                    kernel_l_r = input$l_r,
                    kernel_u_r = input$u_r)) %>% 
                mutate(Numbers = Kernel / 10^((1 + alpha - lambda) * x),
                       Scaled = Numbers / 10^(x * input$s),
                       Biomass = Numbers / 10^(x),
                       Scaled = Scaled / sum(Scaled) / dx,
                       Numbers = Numbers / sum(Numbers) / dx,
                       Biomass = Biomass / sum(Biomass) / dx) %>% 
                gather(key = Type, value = "Density",
                       Numbers, Scaled, Biomass)
            
            pl + geom_line(data = df,
                           aes(x, Density, colour = Type),
                           size = 3) 
        })
        output$plot_kernel <- renderPlot({
            df <- tibble(
                x = x,
                Kernel = power_law_pred_kernel(
                    r, 
                    kernel_exp = input$ex,
                    kernel_l_l = input$l_l,
                    kernel_u_l = input$u_l,
                    kernel_l_r = input$l_r,
                    kernel_u_r = input$u_r)) %>% 
                mutate(Kernel = Kernel / sum(Kernel) / dx)
            
            ggplot(df) + geom_line(aes(x, Kernel)) +
                xlab("Log_10 of predator/prey mass ratio") 
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