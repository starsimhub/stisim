# STIsim HIV Shiny Web App - Simple UI with Tab Panes
# Simplified user interface using tabsetPanel

ui <- fluidPage(
  titlePanel("STIsim HIV Model - Tab Panes Demo"),
  
  # Sidebar with parameters
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Simulation Parameters"),
      
      # Basic simulation parameters
      sliderInput("n_agents", "Population Size", value = 500, min = 100, max = 5000, step = 100),
      sliderInput("start_year", "Start Year", value = 2020, min = 1990, max = 2030, step = 1),
      sliderInput("duration", "Duration (years)", value = 20, min = 1, max = 50, step = 1),
      
      # HIV parameters
      h5("HIV Parameters"),
      sliderInput("init_prev", "Initial HIV Prevalence", value = 0.05, min = 0, max = 1, step = 0.01),
      sliderInput("beta_m2f", "Male-to-Female Transmission Rate", value = 0.05, min = 0, max = 1, step = 0.01),
      sliderInput("beta_m2c", "Mother-to-Child Transmission Rate", value = 0.025, min = 0, max = 1, step = 0.01),
      
      # Action buttons
      br(),
      actionButton("run_simulation", "Run Simulation", class = "btn-primary", width = "100%"),
      br(), br(),
      actionButton("reset_params", "Reset Parameters", class = "btn-secondary", width = "100%")
    ),
    
    # Main content area with tab panels
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        type = "tabs",
        
        # Simulation Tab
        tabPanel("Simulation", icon = icon("play"),
          h3("Simulation Control"),
          p("Configure and run your HIV simulation here."),
          
          # Status box
          wellPanel(
            h4("Simulation Status"),
            textOutput("simulation_status"),
            textOutput("progress_text")
          ),
          
          # Quick results
          fluidRow(
            column(6,
              wellPanel(
                h4("Quick Results"),
                plotlyOutput("quick_plot", height = "300px")
              )
            ),
            column(6,
              wellPanel(
                h4("Population Summary"),
                tableOutput("population_summary")
              )
            )
          )
        ),
        
        # Results Tab
        tabPanel("Results", icon = icon("chart-line"),
          h3("Simulation Results"),
          p("View detailed results and analysis from your simulation."),
          
          # Main prevalence plot
          wellPanel(
            h4("HIV Prevalence Over Time"),
            plotlyOutput("prevalence_plot", height = "400px")
          ),
          
          # Secondary plots
          fluidRow(
            column(6,
              wellPanel(
                h4("HIV Incidence"),
                plotlyOutput("incidence_plot", height = "300px")
              )
            ),
            column(6,
              wellPanel(
                h4("CD4 Count Distribution"),
                plotlyOutput("cd4_plot", height = "300px")
              )
            )
          ),
          
          fluidRow(
            column(6,
              wellPanel(
                h4("Treatment Coverage"),
                plotlyOutput("treatment_plot", height = "300px")
              )
            ),
            column(6,
              wellPanel(
                h4("Age-specific Prevalence"),
                plotlyOutput("age_prev_plot", height = "300px")
              )
            )
          ),
          
          # Network analysis
          wellPanel(
            h4("Network Analysis"),
            plotlyOutput("network_plot", height = "300px")
          )
        ),
        
        # Parameters Tab
        tabPanel("Parameters", icon = icon("cog"),
          h3("Parameter Analysis"),
          p("Review current parameters and sensitivity analysis."),
          
          wellPanel(
            h4("Current Parameters"),
            DTOutput("parameters_table")
          ),
          
          wellPanel(
            h4("Parameter Sensitivity"),
            plotlyOutput("sensitivity_plot", height = "400px")
          )
        ),
        
        # About Tab
        tabPanel("About", icon = icon("info"),
          h3("About STIsim HIV Model"),
          
          wellPanel(
            h4("STIsim HIV Web Application"),
            p("This web application provides an interactive interface for the STIsim HIV modeling framework."),
            p("STIsim is an agent-based modeling framework for simulating sexually transmitted diseases, built on the Starsim architecture."),
            
            h5("Features:"),
            tags$ul(
              tags$li("Agent-based HIV simulation"),
              tags$li("Configurable disease parameters"),
              tags$li("Sexual network modeling"),
              tags$li("Intervention modeling (testing, treatment, prevention)"),
              tags$li("Interactive visualizations"),
              tags$li("Parameter sensitivity analysis")
            ),
            
            h5("Model Components:"),
            tags$ul(
              tags$li("HIV natural history (acute, latent, late-stage)"),
              tags$li("Sexual transmission networks"),
              tags$li("Demographic processes"),
              tags$li("HIV testing and treatment"),
              tags$li("Prevention interventions")
            ),
            
            h5("Contact:"),
            p("For questions or support, please visit the STIsim GitHub repository or contact info@starsim.org"),
            
            h5("Version:"),
            p("STIsim v1.4 - Shiny Web Interface v1.0")
          )
        )
      )
    )
  )
)