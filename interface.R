#interface.R                   

#-------------------------------------------------------------------------------
#
#  Script for generating the Interface and overall appearance
#
# Copyright (C) 2025 Jean-Pierre Gnimatin, Marlon Grodd, Susanne Weber, Derek Hazard, Martin Wolkewitz
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------


# Defining the UI of the RShiny App------------------------------------------------------------------------------------------------------


ui <- fluidPage(
  tags$i("Copyright © 2025 Jean-Pierre Gnimatin, Marlon Grodd, Susanne Weber, Derek Hazard, Martin Wolkewitz"),
  
  # Adding logos at the top, centered
  tags$div(
    style = "text-align: center; margin-bottom: 0px;",
    tags$img(src = "fdm_logo.png", height = "150px", style = "margin-right: 300px;"),
    tags$img(src = "ukl_logo.png", height = "90px")
  ),
  
  
  
  # Shiny App Tilte
  titlePanel("HAISim (Hospital-Acquired infection Interventions Simulator) Shiny App"),
  
  tags$p("This Shiny App simulates and evaluates the impact of various intervention strategies on HAIs outcomes by modeling different intervention factors.", 
         style = "text-align: left; font-weight: bold;"),
  
  tags$p("States: 0- Admission into hospital, 1- Hospital-Acquired Infection, 2-Discharged alive, 3-Death, 4-Discharged Alive after HAI, 5-Death after HAI", 
         style = "text-align: left; font-weight: bold;"),
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("setting", "Select Setting:", choices = c("Setting 1 (Enhanced Treatment)", "Setting 2 (Enhanced Treatment + Prevention)")),
      
      # Setting 1 inputs
      conditionalPanel(
        condition = "input.setting == 'Setting 1 (Enhanced Treatment)'",
        h4("Setting 1 (Enhanced Treatment): Parameter Inputs"),
        numericInput("lambda01_m1", "λ01 (Infection hazard rate):", value = 0.01, min = 0.001, step = 0.001),
        numericInput("lambda02_m1", "λ02 (Discharge hazard rate without infection):", value = 0.2, min = 0.001, step = 0.001),
        numericInput("lambda03_m1", "λ03 (Death hazard rate without infection):", value = 0.05, min = 0.001, step = 0.001),
        numericInput("lambda14_m1", "λ14 (Discharge hazard rate with infection):", value = 0.1, min = 0.001, step = 0.001),
        numericInput("lambda15_m1", "λ15 (Death hazard rate with infection):", value = 0.05, min = 0.001, step = 0.001),
        sliderInput("alpha_range_m1", "Alpha Range (Impacting Discharge after HAI):", min = 0, max = 1, value = c(0, 1), step = 0.02),
        sliderInput("beta_range_m1", "Beta Range  (Impacting Death after HAI):", min = 0, max = 1, value = c(0, 1), step = 0.02),
        numericInput("n_patients_m1", "Total Number of Patients (n):", value = 10000, min = 1, step = 1),
        actionButton("simulate_m1", "Run Setting 1 Simulation")
      ),
      
      # Setting 2 inputs
      conditionalPanel(
        condition = "input.setting == 'Setting 2 (Enhanced Treatment + Prevention)'",
        h4("Setting 2 (Enhanced Treatment + Prevention): Lambda Parameter Inputs"),
        numericInput("lambda01_m2", "λ01 (Infection hazard rate):", value = 0.01, min = 0),
        numericInput("lambda02_m2", "λ02 (Discharge hazard rate without infection):", value = 0.2, min = 0),
        numericInput("lambda03_m2", "λ03 (Death hazard rate without infection):", value = 0.05, min = 0),
        numericInput("lambda14_m2", "λ14 (Discharge hazard rate with infection):", value = 0.2, min = 0),
        numericInput("lambda15_m2", "λ15 (Death hazard rate with infection):", value = 0.25, min = 0),
        h4("Simulation Grid"),
        sliderInput("alpha_range_m2", "Alpha Range (Impacting Discharge after HAI):", min = 0.1, max = 1, value = c(0.1, 1), step = 0.1),
        sliderInput("beta_range_m2", "Beta Range (Impacting Death after HAI):", min = 0.1, max = 1, value = c(0.1, 1), step = 0.1),
        sliderInput("Theta_range_m2", "Theta Range (Impacting Prevention of HAI):", min = 0.1, max = 1, value = c(0.1, 1), step = 0.1),
        numericInput("num_patients_m2", "Total Number of Patients (n):", value = 10000, min = 1, step = 1),
        actionButton("simulate_m2", "Run Setting 2 Simulation")
      )
    ),
    
    # Main panel
    mainPanel(
      conditionalPanel(
        condition = "input.setting == 'Setting 1 (Enhanced Treatment)'",
        fluidRow(
          column(6,
                 tags$div(
                   style = "padding: 10px;",
                   h4("About Setting 1"),
                   p("Setting 1 evaluates the impact of ENHANCED TREATMENT interventions (against HAI) on patient outcomes across varying α and β values."),
                   p("α is the effect of ENHANCED TREATMENT interventions on DISCHARGE and β on DEATH (See formulas attached to Setting 1 Figure)."),
                   p("α and β range from 0 to 1, where 1 indicates no treatment effect (minimal effect) and 0 represents the maximum treatment effect."),
                   p("Transition Hazard Rate (λᵢⱼ): This defines the rate at which patients move from state i to state j."),
                   p("λᵢⱼ = (Number of transitions from state i to state j) ÷ (Total patient-days in state i)."),
                 )
          ),
          column(6, 
                 tags$div(
                   style = "text-align: center;",
                   tags$img(src = "Setting1.jpg", height = "170px")
                 )
          )
          
        ),
        tabsetPanel(
          tabPanel("Reduced Mortality", "Reduced mortality rate calculated for specified α and β values and displayed in the plots (e.g. RM=0.05)", plotlyOutput("RM_plot_m1", height = "700px")),
          
          
          tabPanel("Lives Saved", "Number of lives saved per population unit, calculated for specified α and β values and displayed in the plots (e.g. 50 per 1,000 patients)  / Note: Lives Saved = Reduced Mortality × Total number of Patients",plotlyOutput("Lives_Saved_plot_m1", height = "700px")),
          
          
          tabPanel("Change in Length of Stay", "Change in hospital length of stay relative to baseline, calculated for specified α and β values and displayed in the plots (e.g. C_LoS=2.5days)", plotlyOutput("C_LoS_plot_m1", height = "700px")),
          
          
          tabPanel("Patient-Days Reduced", "Total reduction in patient‑days calculated for specified α and β values and displayed in the plots (e.g. 250 patient‑days) / Note: Patient-Days Reduced = Change in Length of Stay × Total number of Patients", plotlyOutput("Patient_Days_Reduced_plot_m1", height = "700px"))
        )),
      
      
      
      
      
      conditionalPanel(
        condition = "input.setting == 'Setting 2 (Enhanced Treatment + Prevention)'",
        fluidRow(
          column(6,
                 tags$div(
                   style = "padding: 20px;",            
                   h4("About Setting 2"),
                   p("Setting 2 evaluates the impact of ENHANCED TREATMENT and INFECTION PREVENTION interventions on patient outcomes across varying α, β and θ values."),
                   p("α is the effect of ENHANCED TREATMENT interventions on DISCHARGE and β on DEATH."),
                   p("θ is the IMPROVED INFECTION PREVENTION factor."),
                   p("α, β and θ range from 0 to 1, where 1 indicates no treatment effect (minimal effect) and 0 represents the maximum treatment effect."),
                   p("Transition Hazard Rate (λᵢⱼ): This defines the rate at which patients move from state i to state j."),
                   p("λᵢⱼ = (Number of transitions from state i to state j) ÷ (Total patient-days in state i)."),
                 )
          ),
          column(6, 
                 tags$div(
                   style = "text-align: center;",
                   tags$img(src = "Setting2.jpg", height = "170px")
                 )
          )
          
        ),
        tabsetPanel(
          tabPanel("Reduced Mortality", plotlyOutput("RM_plot_m2", height = "650px")),
          tabPanel("Lives Saved", plotlyOutput("Lives_Saved_plot_m2", height = "650px")),
          tabPanel("Change in Length of Stay", plotlyOutput("C_LoS_plot_m2", height = "650px")),
          tabPanel("Patient-Days Reduced", plotlyOutput("Patient_Days_Reduced_plot_m2", height = "650px"))
        )
      )
    )
  )
)






# Define Server --------------------------------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  ### Server Setting 1
  
  observeEvent(input$simulate_m1, {
    alpha_values_m1 <- seq(input$alpha_range_m1[1], input$alpha_range_m1[2], by = 0.02)
    beta_values_m1 <- seq(input$beta_range_m1[1], input$beta_range_m1[2], by = 0.02)
    grid_m1 <- expand.grid(alpha_m1 = alpha_values_m1, beta_m1 = beta_values_m1)
    
    
    
    # Calculation Function
    Calculations_m1 <- function(alpha_m1, beta_m1, lambda02_m1, lambda03_m1, lambda14_m1, lambda15_m1, lambda01_m1) {
      lambda14T_m1 <- alpha_m1 * (lambda14_m1 - lambda02_m1) + lambda02_m1
      lambda15T_m1 <- beta_m1 * (lambda15_m1 - lambda03_m1) + lambda03_m1
      lambda1T_m1 <- lambda14T_m1 + lambda15T_m1
      
      lambda0_m1 <- lambda01_m1 + lambda02_m1 + lambda03_m1
      lambda1_m1 <- lambda14_m1 + lambda15_m1
      
      Prob_0_m1 <- (lambda1_m1 * lambda03_m1 + lambda01_m1 * lambda15_m1) / (lambda0_m1 * lambda1_m1)
      Prob_1_m1 <- (lambda1T_m1 * lambda03_m1 + lambda01_m1 * lambda15T_m1) / (lambda0_m1 * lambda1T_m1)
      RM_m1 <- Prob_0_m1 - Prob_1_m1  
      
      LoS_0_m1 <- (lambda1_m1 + lambda01_m1) / (lambda0_m1 * lambda1_m1)
      LoS_1_m1 <- (lambda1T_m1 + lambda01_m1) / (lambda0_m1 * lambda1T_m1)
      C_LoS_m1 <- LoS_0_m1 - LoS_1_m1
      
      return(c(RM_m1, C_LoS_m1))
    }
    
    
    
    
    # Performing Calculations
    results_m1 <- t(mapply(Calculations_m1, grid_m1$alpha_m1, grid_m1$beta_m1, 
                           MoreArgs = list(lambda02_m1 = input$lambda02_m1,
                                           lambda03_m1 = input$lambda03_m1,
                                           lambda14_m1 = input$lambda14_m1,
                                           lambda15_m1 = input$lambda15_m1,
                                           lambda01_m1 = input$lambda01_m1)))
    
    grid_m1$RM_m1 <- results_m1[, 1]  
    grid_m1$C_LoS_m1 <- results_m1[, 2]
    
    n_m1 <- input$n_patients_m1
    grid_m1$Lives_Saved_m1 <- grid_m1$RM_m1 * n_m1
    grid_m1$Patient_Days_Reduced_m1 <- grid_m1$C_LoS_m1 * n_m1 
    
    
    # Plot Outputs for Setting 1
    
    #RM
    
    grid_m1 <- grid_m1 %>% rename(RM = RM_m1)
    
    
    output$RM_plot_m1 <- renderPlotly({
      plot_ly(grid_m1, x = ~alpha_m1, y = ~beta_m1, z = ~RM, 
              color = ~RM, colors = viridis(100),
              type = "scatter3d", mode = "markers", 
              marker = list(size = 3),
              text = ~paste0("Alpha: ", round(alpha_m1, 2), 
                             "<br>Beta: ", round(beta_m1, 2), 
                             "<br>RM: ", round(RM, 4)),
              hoverinfo = "text",
              name= "Reduced M") %>%
        layout(
          scene = list(xaxis = list(title = "Alpha"),
                       yaxis = list(title = "Beta"),
                       zaxis = list(title = "Reduced Mortality (RM)"),
                       camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
                       aspectratio = list(x = 1.5, y = 1.5, z = 1)))
    })
    
    
    
    #Change in Length of Stay (C_LOS)
    
    grid_m1 <- grid_m1 %>% rename(C_LoS = C_LoS_m1)
    
    output$C_LoS_plot_m1 <- renderPlotly({
      plot_ly(grid_m1, x = ~alpha_m1, y = ~beta_m1, z = ~C_LoS, 
              color = ~C_LoS, colors = viridis(100),
              type = "scatter3d", mode = "markers", 
              marker = list(size = 3),
              text = ~paste0("Alpha: ", round(alpha_m1, 2), 
                             "<br>Beta: ", round(beta_m1, 2), 
                             "<br>C_LoS: ", round(C_LoS, 4)),
              hoverinfo = "text") %>%
        layout(
          scene = list(xaxis = list(title = "Alpha"),
                       yaxis = list(title = "Beta"),
                       zaxis = list(title = "Change in Length of Stay"),
                       camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
                       aspectratio = list(x = 1.5, y = 1.5, z = 1)))
    })
    
    
    #Lives Saved
    
    grid_m1 <- grid_m1 %>% rename(Lives_Saved = Lives_Saved_m1)
    
    output$Lives_Saved_plot_m1 <- renderPlotly({
      plot_ly(grid_m1, x = ~alpha_m1, y = ~beta_m1, z = ~Lives_Saved, 
              color = ~Lives_Saved, colors = viridis(100),
              type = "scatter3d", mode = "markers", 
              marker = list(size = 3),
              text = ~paste0("Alpha: ", round(alpha_m1, 2), 
                             "<br>Beta: ", round(beta_m1, 2), 
                             "<br>Lives Saved: ", round(Lives_Saved, 4)),
              hoverinfo = "text") %>%
        layout(
          scene = list(xaxis = list(title = "Alpha"),
                       yaxis = list(title = "Beta"),
                       zaxis = list(title = "Lives Saved"),
                       camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
                       aspectratio = list(x = 1.5, y = 1.5, z = 1)))
    })
    
    # Reduced Patient-days
    
    grid_m1 <- grid_m1 %>% rename(Patient_Days_Reduced = Patient_Days_Reduced_m1)
    
    
    output$Patient_Days_Reduced_plot_m1 <- renderPlotly({
      plot_ly(grid_m1, x = ~alpha_m1, y = ~beta_m1, z = ~Patient_Days_Reduced, 
              color = ~Patient_Days_Reduced, colors = viridis(100),
              type = "scatter3d", mode = "markers", 
              marker = list(size = 3),
              text = ~paste0("Alpha: ", round(alpha_m1, 2), 
                             "<br>Beta: ", round(beta_m1, 2), 
                             "<br>Patient-Days Reduced: ", round(Patient_Days_Reduced, 4)),
              hoverinfo = "text") %>%
        layout(
          scene = list(xaxis = list(title = "Alpha"),
                       yaxis = list(title = "Beta"),
                       zaxis = list(title = "Patient-Days Reduced"),
                       camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
                       aspectratio = list(x = 1.5, y = 1.5, z = 1)))
    })
  })
  
  
  
  
  
  
  
  
  
  
  ### Server Setting 2
  
  simulation_results_m2 <- eventReactive(input$simulate_m2, {
    
    
    
    # Retrieve lambda values from user input
    lambda_values_m2 <- list(
      lambda01_m2 = input$lambda01_m2,
      lambda02_m2 = input$lambda02_m2,
      lambda03_m2 = input$lambda03_m2,
      lambda14_m2 = input$lambda14_m2,
      lambda15_m2 = input$lambda15_m2
    )
    
    
    
    
    # Generating sequences for alpha, beta, and Theta
    alpha_values_m2 <- seq(input$alpha_range_m2[1], input$alpha_range_m2[2], length.out = 10)
    beta_values_m2 <- seq(input$beta_range_m2[1], input$beta_range_m2[2], length.out = 10)
    Theta_values_m2 <- seq(input$Theta_range_m2[1], input$Theta_range_m2[2], length.out = 10)
    
    
    
    
    # Creating a grid of all combinations
    grid_m2 <- expand.grid(alpha_m2 = alpha_values_m2, beta_m2 = beta_values_m2, Theta = Theta_values_m2)
    
    
    
    # Function to compute metrics
    calculate_metrics_m2 <- function(alpha_m2, beta_m2, Theta, lambda_values_m2) {
      
      # Adjust lambda values based on alpha and beta
      lambda14T_m2 <- alpha_m2 * (lambda_values_m2$lambda14_m2 - lambda_values_m2$lambda02_m2) + lambda_values_m2$lambda02_m2
      lambda15T_m2 <- beta_m2 * (lambda_values_m2$lambda15_m2 - lambda_values_m2$lambda03_m2) + lambda_values_m2$lambda03_m2
      
      lambda1T_m2 <- lambda14T_m2 + lambda15T_m2
      lambda0_m2 <- lambda_values_m2$lambda01_m2 + lambda_values_m2$lambda02_m2 + lambda_values_m2$lambda03_m2
      lambda1_m2 <- lambda_values_m2$lambda14_m2 + lambda_values_m2$lambda15_m2
      
      
      
      # Calculate probabilities and metrics
      Prob_0_m2 <- (lambda1_m2 * lambda_values_m2$lambda03_m2 + lambda_values_m2$lambda01_m2 * lambda_values_m2$lambda15_m2) / (lambda0_m2 * lambda1_m2)
      Prob_2_m2 <- (lambda1T_m2 * lambda_values_m2$lambda03_m2 + Theta * lambda_values_m2$lambda01_m2 * lambda15T_m2) / (lambda0_m2 * lambda1T_m2)
      
      RM_m2 <- Prob_0_m2 - Prob_2_m2
      LoS_0_m2 <- (lambda1_m2 + lambda_values_m2$lambda01_m2) / (lambda0_m2 * lambda1_m2)
      LoS_2_m2 <- (Theta * lambda_values_m2$lambda01_m2 + lambda1T_m2) / (lambda0_m2 * lambda1T_m2)
      C_LoS_m2 <- LoS_0_m2 - LoS_2_m2
      
      return(list(RM_m2 = RM_m2, C_LoS_m2 = C_LoS_m2))
    }
    
    
    
    # Applying the function to the grid and storing results
    results <- mapply(calculate_metrics_m2, grid_m2$alpha_m2, grid_m2$beta_m2, grid_m2$Theta,
                      MoreArgs = list(lambda_values_m2 = lambda_values_m2))
    
    
    
    # Binding results to the grid
    grid_m2$RM_m2 <- unlist(results[1, ])
    grid_m2$C_LoS_m2 <- unlist(results[2, ])
    grid_m2$Lives_Saved_m2 <- grid_m2$RM_m2 * input$num_patients_m2
    grid_m2$Patient_Days_Reduced_m2 <- grid_m2$C_LoS_m2 * input$num_patients_m2
    
    
    return(grid_m2)
  })
  
  
  
  
  
  # Plots for Setting 2 (Optimized with native plotly) ----------------------------
  # Load necessary libraries
  library(plotly)
  library(viridis)
  
  # Modified generate_plots function with theta rounding
  generate_plots <- function(data, z_col, plot_title, colorbar_title, hover_template) {
    
    # Add rounding to handle floating-point precision
    data$Theta_rounded <- round(data$Theta, 1)  # Round to 1 decimal place
    
    # Split data by rounded theta values
    theta_groups <- split(data, data$Theta_rounded)
    
    # Determine number of plots
    num_plots <- length(theta_groups)
    
    # Determine rows
    n_rows <- ceiling(num_plots / 5)
    
    # Create list of subplots
    plots <- lapply(seq_along(theta_groups), function(i) {
      theta_group <- theta_groups[[i]]
      
      # Check if it's in the 2nd row
      is_second_row <- (i > 5)
      
      plot_ly(
        data = theta_group,
        x = ~alpha_m2,
        y = ~beta_m2,
        z = as.formula(paste0("~", z_col)),
        type = "heatmap",
        colors = viridis_pal(option = "viridis")(256),
        zmin = min(data[[z_col]], na.rm = TRUE),
        zmax = max(data[[z_col]], na.rm = TRUE),
        zauto = FALSE,
        hoverinfo = "x+y+z",
        hovertemplate = hover_template,
        showscale = FALSE,
        colorbar = list(title = "")  
      ) %>%
        layout(
          xaxis = list(
            title = "Alpha",
            showgrid = FALSE,
            tickformat = ".2f",
            tickvals = if (is_second_row) c(0.2, 0.4, 0.6, 0.8, 1) else NULL
          ),
          yaxis = list(
            title = "Beta",
            showgrid = FALSE,
            tickformat = ".2f"
          ),
          annotations = list(
            text = paste("θ =", unique(theta_group$Theta_rounded)),  # Use rounded value
            font = list(size = 12),
            x = 0.5,
            y = if (is_second_row) 1.02 else 1.06,
            xref = "paper",
            yref = "paper",
            showarrow = FALSE
          ),
          margin = list(t = 50, b = 40)
        )
    })
    
    # Create shared colorbar
    colorbar <- list(
      title = NULL,
      thickness = 20,
      x = 1.02,
      len = 0.8
    )
    
    # Combine subplots
    p <- subplot(
      plots,
      nrows = n_rows,
      margin = 0.015,
      shareX = TRUE,
      shareY = TRUE
    ) %>%
      layout(
        title = list(
          text = plot_title,
          y = 0.99,
          font = list(size = 18)
        ),
        coloraxis = list(colorbar = colorbar),
        margin = list(autoexpand = TRUE, pad = 5),
        autosize = TRUE
      )
    
    # Add colorbar to first plot only
    p$x$data[[1]]$showscale <- TRUE
    p
  }
  
  # --------------------- Plot Outputs -----------------------------
  
  # Reduced Mortality (RM2) Plot
  output$RM_plot_m2 <- renderPlotly({
    data <- simulation_results_m2()
    generate_plots(
      data,
      z_col = "RM_m2",
      plot_title = "Reduced Mortality (RM)",
      colorbar_title = "RM",
      hover_template = paste(
        "Alpha: %{x:.2f}<br>",
        "Beta: %{y:.2f}<br>",
        "RM: %{z:.4f}<extra></extra>"
      )
    )
  })
  
  # Lives Saved Plot2
  output$Lives_Saved_plot_m2 <- renderPlotly({
    data <- simulation_results_m2()
    generate_plots(
      data,
      z_col = "Lives_Saved_m2",
      plot_title = "Lives Saved (RM × Total number of Patients)",
      colorbar_title = "Lives Saved",
      hover_template = paste(
        "Alpha: %{x:.2f}<br>",
        "Beta: %{y:.2f}<br>",
        "Lives Saved: %{z:.1f}<extra></extra>"
      )
    )
  })
  
  # Change in Length of Stay (C_LoS2) Plot
  output$C_LoS_plot_m2 <- renderPlotly({
    data <- simulation_results_m2()
    generate_plots(
      data,
      z_col = "C_LoS_m2",
      plot_title = "Change in Length of Stay (C_LoS)",
      colorbar_title = "C_LoS",
      hover_template = paste(
        "Alpha: %{x:.2f}<br>",
        "Beta: %{y:.2f}<br>",
        "C_LoS: %{z:.4f}<extra></extra>"
      )
    )
  })
  
  # Patient-Days Reduced Plot2
  output$Patient_Days_Reduced_plot_m2 <- renderPlotly({
    data <- simulation_results_m2()
    generate_plots(
      data,
      z_col = "Patient_Days_Reduced_m2",
      plot_title = "Patient-Days Reduced (C_LoS × Total number of Patients)",
      colorbar_title = "Patient-Days",
      hover_template = paste(
        "Alpha: %{x:.2f}<br>",
        "Beta: %{y:.2f}<br>",
        "Days Reduced: %{z:.1f}<extra></extra>"
      )
    )
  })
  
  
  
}


# Run Shiny App ------------------------------------------------------------------------------------------------------------

shinyApp(ui, server)