# STIsim HIV Shiny Web App - Server
# Server logic and reactive functions

server <- function(input, output, session) {
  
  # Generate mock simulation data when Python is not available
  generate_mock_data <- function(params) {
    years <- params$start:(params$start + params$duration)
    n_years <- length(years)
    
    # Generate realistic HIV prevalence curve
    base_prev <- params$init_prev
    growth_rate <- params$beta_m2f * 0.1  # Scale down for realistic growth
    prevalence <- base_prev * (1 + growth_rate * (0:(n_years-1)))
    prevalence <- pmin(prevalence, 0.3)  # Cap at 30%
    
    # Generate incidence (new infections)
    incidence <- diff(c(0, prevalence)) * params$n_agents
    incidence[incidence < 0] <- 0
    
    # Generate CD4 counts for infected individuals
    n_infected <- round(prevalence[n_years] * params$n_agents)
    cd4_counts <- if (n_infected > 0) {
      rnorm(n_infected, mean = 400, sd = 150)
    } else {
      numeric(0)
    }
    
    # Generate age-specific prevalence
    age_groups <- c('15-24', '25-34', '35-44', '45-54', '55+')
    age_prev <- setNames(rep(0, length(age_groups)), age_groups)
    if (n_infected > 0) {
      age_prev <- setNames(runif(length(age_groups), 0, prevalence[n_years] * 2), age_groups)
    }
    
    # Network statistics
    network_stats <- list(
      mean_partners = runif(1, 1, 3),
      max_partners = sample(5:20, 1)
    )
    
    return(list(
      timevec = years,
      years = years,
      hiv_prevalence = prevalence,
      hiv_incidence = incidence,
      population_summary = list(
        total_pop = params$n_agents,
        hiv_infected = n_infected,
        hiv_prevalence = prevalence[n_years],
        on_art = round(n_infected * 0.6),  # 60% on ART
        art_coverage = 0.6
      ),
      age_prevalence = age_prev,
      cd4_counts = cd4_counts,
      network_stats = network_stats,
      parameters = params
    ))
  }
  
  # Reactive values to store simulation results
  simulation_results <- reactiveValues(
    data = NULL,
    status = "Ready to run simulation",
    progress = 0
  )
  
  # Parameter validation
  validate_parameters <- function() {
    errors <- c()
    
    if (input$n_agents < 100) {
      errors <- c(errors, "Population size must be at least 100")
    }
    
    if (input$init_prev < 0 || input$init_prev > 1) {
      errors <- c(errors, "Initial prevalence must be between 0 and 1")
    }
    
    if (input$beta_m2f < 0 || input$beta_m2f > 1) {
      errors <- c(errors, "Male-to-female transmission rate must be between 0 and 1")
    }
    
    if (input$duration < 1) {
      errors <- c(errors, "Duration must be at least 1 year")
    }
    
    return(errors)
  }
  
  # Run simulation
  observeEvent(input$run_simulation, {
    # Validate parameters
    errors <- validate_parameters()
    if (length(errors) > 0) {
      showNotification(paste("Parameter errors:", paste(errors, collapse = ", ")), 
                     type = "error", duration = 5)
      return()
    }
    
    # Update status
    simulation_results$status <- "Running simulation..."
    simulation_results$progress <- 10
    
    # Prepare parameters
    params <- list(
      n_agents = as.integer(input$n_agents),
      start = input$start_year,
      dur = input$duration,
      init_prev = input$init_prev,
      beta_m2f = input$beta_m2f,
      beta_m2c = input$beta_m2c,
      beta_m2m = input$beta_m2m,
      rel_beta_f2m = input$rel_beta_f2m,
      eff_condom = input$eff_condom,
      dur_acute = input$dur_acute,
      dur_latent = input$dur_latent,
      dur_falling = input$dur_falling,
      rel_trans_acute = input$rel_trans_acute,
      rel_trans_falling = input$rel_trans_falling,
      include_testing = input$include_testing,
      include_art = input$include_art,
      include_vmmc = input$include_vmmc,
      include_prep = input$include_prep,
      network_type = input$network_type,
      debut_age = input$debut_age
    )
    
    # Run simulation with progress updates
    tryCatch({
      simulation_results$progress <- 30
      
      # Call Python simulation function or use mock data
      if (exists("python_available") && python_available) {
        result <- py$run_hiv_simulation(params)
      } else {
        # Generate mock data when Python is not available
        result <- generate_mock_data(params)
      }
      
      simulation_results$progress <- 70
      
      # Store results
      simulation_results$data <- result
      if (exists("python_available") && python_available) {
        simulation_results$status <- "Simulation completed successfully"
      } else {
        simulation_results$status <- "Simulation completed (mock data - Python not available)"
      }
      simulation_results$progress <- 100
      
      showNotification("Simulation completed successfully!", type = "message")
      
    }, error = function(e) {
      simulation_results$status <- paste("Simulation failed:", e$message)
      simulation_results$progress <- 0
      showNotification(paste("Simulation failed:", e$message), type = "error")
    })
  })
  
  # Reset parameters
  observeEvent(input$reset_params, {
    updateSliderInput(session, "n_agents", value = 500)
    updateSliderInput(session, "start_year", value = 2020)
    updateSliderInput(session, "duration", value = 20)
    updateSliderInput(session, "init_prev", value = 0.05)
    updateSliderInput(session, "beta_m2f", value = 0.05)
    updateSliderInput(session, "beta_m2c", value = 0.025)
    updateSliderInput(session, "beta_m2m", value = 0.1)
    updateSliderInput(session, "rel_beta_f2m", value = 0.5)
    updateSliderInput(session, "eff_condom", value = 0.9)
    updateSliderInput(session, "dur_acute", value = 3)
    updateSliderInput(session, "dur_latent", value = 10)
    updateSliderInput(session, "dur_falling", value = 3)
    updateSliderInput(session, "rel_trans_acute", value = 6)
    updateSliderInput(session, "rel_trans_falling", value = 8)
    updateCheckboxInput(session, "include_testing", value = TRUE)
    updateCheckboxInput(session, "include_art", value = TRUE)
    updateCheckboxInput(session, "include_vmmc", value = TRUE)
    updateCheckboxInput(session, "include_prep", value = FALSE)
    updateSelectInput(session, "network_type", selected = "structured")
    updateSliderInput(session, "debut_age", value = 20)
    
    simulation_results$data <- NULL
    simulation_results$status <- "Parameters reset"
    simulation_results$progress <- 0
  })
  
  # Simulation status output
  output$simulation_status <- renderText({
    simulation_results$status
  })
  
  # Progress bar - removed updateProgressBar call as it's not available
  # Progress is shown via text output instead
  
  # Population summary
  output$population_summary <- renderTable({
    if (!is.null(simulation_results$data)) {
      data.frame(
        Metric = c("Total Population", "HIV Infected", "HIV Prevalence", "On ART", "ART Coverage"),
        Value = c(
          simulation_results$data$population_summary$total_pop,
          simulation_results$data$population_summary$hiv_infected,
          paste0(round(simulation_results$data$population_summary$hiv_prevalence * 100, 2), "%"),
          simulation_results$data$population_summary$on_art,
          paste0(round(simulation_results$data$population_summary$art_coverage * 100, 2), "%")
        )
      )
    } else {
      data.frame(Metric = "No simulation data", Value = "Run simulation first")
    }
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Quick plot
  output$quick_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      plot_ly(x = data$years, y = data$hiv_prevalence * 100, type = 'scatter', mode = 'lines+markers',
              line = list(color = 'red', width = 3),
              marker = list(color = 'red', size = 4)) %>%
        layout(title = 'HIV Prevalence Over Time',
               xaxis = list(title = 'Year'),
               yaxis = list(title = 'Prevalence (%)'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Prevalence plot
  output$prevalence_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      plot_ly(x = data$years, y = data$hiv_prevalence * 100, type = 'scatter', mode = 'lines',
              line = list(color = 'red', width = 3)) %>%
        layout(title = 'HIV Prevalence Over Time',
               xaxis = list(title = 'Year'),
               yaxis = list(title = 'Prevalence (%)'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Incidence plot
  output$incidence_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      plot_ly(x = data$years, y = data$hiv_incidence, type = 'scatter', mode = 'lines',
              line = list(color = 'orange', width = 3)) %>%
        layout(title = 'HIV Incidence Over Time',
               xaxis = list(title = 'Year'),
               yaxis = list(title = 'New Infections'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # CD4 plot
  output$cd4_plot <- renderPlotly({
    if (!is.null(simulation_results$data) && length(simulation_results$data$cd4_counts) > 0) {
      data <- simulation_results$data
      plot_ly(x = data$cd4_counts, type = 'histogram', 
              marker = list(color = 'lightblue')) %>%
        layout(title = 'CD4 Count Distribution (HIV Infected)',
               xaxis = list(title = 'CD4 Count (cells/μL)'),
               yaxis = list(title = 'Frequency'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No CD4 data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Treatment plot
  output$treatment_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      art_coverage <- seq(0, data$population_summary$art_coverage, length.out = length(data$years))
      
      plot_ly(x = data$years, y = art_coverage * 100, type = 'scatter', mode = 'lines',
              line = list(color = 'green', width = 3)) %>%
        layout(title = 'ART Coverage Over Time',
               xaxis = list(title = 'Year'),
               yaxis = list(title = 'Coverage (%)'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Age-specific prevalence plot
  output$age_prev_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      age_groups <- names(data$age_prevalence)
      prevalences <- unlist(data$age_prevalence) * 100
      
      plot_ly(x = age_groups, y = prevalences, type = 'bar',
              marker = list(color = 'purple')) %>%
        layout(title = 'HIV Prevalence by Age Group',
               xaxis = list(title = 'Age Group'),
               yaxis = list(title = 'Prevalence (%)'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Network plot
  output$network_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      network_stats <- data$network_stats
      
      plot_ly(x = c('Mean Partners', 'Max Partners'), 
              y = c(network_stats$mean_partners, network_stats$max_partners),
              type = 'bar', marker = list(color = 'teal')) %>%
        layout(title = 'Sexual Network Statistics',
               xaxis = list(title = 'Metric'),
               yaxis = list(title = 'Number of Partners'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
  
  # Parameters table
  output$parameters_table <- renderDT({
    if (!is.null(simulation_results$data)) {
      params_df <- data.frame(
        Parameter = names(simulation_results$data$parameters),
        Value = unlist(simulation_results$data$parameters)
      )
      datatable(params_df, options = list(pageLength = 20, scrollX = TRUE))
    } else {
      data.frame(Parameter = "No simulation data", Value = "Run simulation first")
    }
  })
  
  # Sensitivity plot
  output$sensitivity_plot <- renderPlotly({
    if (!is.null(simulation_results$data)) {
      data <- simulation_results$data
      params <- data$parameters
      param_names <- names(params)
      param_values <- unlist(params)
      
      # Normalize values for display
      normalized_values <- sapply(param_values, function(val) {
        if (is.numeric(val)) {
          if (val > 1) val / 10 else val * 10
        } else 1
      })
      
      plot_ly(x = param_names, y = normalized_values, type = 'bar',
              marker = list(color = 'lightcoral')) %>%
        layout(title = 'Parameter Values (Normalized)',
               xaxis = list(title = 'Parameters', tickangle = -45),
               yaxis = list(title = 'Normalized Values'))
    } else {
      plot_ly() %>% 
        add_annotations(text = "No simulation data available", 
                       x = 0.5, y = 0.5, showarrow = FALSE) %>%
        layout(xaxis = list(showticklabels = FALSE, showgrid = FALSE),
               yaxis = list(showticklabels = FALSE, showgrid = FALSE))
    }
  })
}
