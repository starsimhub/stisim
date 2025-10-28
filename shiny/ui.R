# STIsim HIV Shiny Web App - UI
# User interface definition

ui <- dashboardPage(
  dashboardHeader(title = "STIsim HIV Model"),
  
  dashboardSidebar(
    width = 350,
    
    # Parameter controls
    div(id = "parameter_controls",
        h4("Simulation Parameters"),
        
        # Basic simulation parameters
        div(class = "slider-container",
            sliderInput("n_agents", "Population Size", value = 500, min = 100, max = 50000, step = 100, 
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("start_year", "Start Year", value = 2020, min = 1990, max = 2030, step = 1,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("duration", "Duration (years)", value = 20, min = 1, max = 50, step = 1,
                       ticks = TRUE, width = "100%")
        ),
        
        # HIV disease parameters
        h5("HIV Parameters"),
        div(class = "slider-container",
            sliderInput("init_prev", "Initial HIV Prevalence", value = 0.05, min = 0, max = 1, step = 0.01,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("beta_m2f", "Male-to-Female Transmission Rate", value = 0.05, min = 0, max = 1, step = 0.01,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("beta_m2c", "Mother-to-Child Transmission Rate", value = 0.025, min = 0, max = 1, step = 0.01,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("beta_m2m", "Male-to-Male Transmission Rate", value = 0.1, min = 0, max = 1, step = 0.01,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("rel_beta_f2m", "Female-to-Male Relative Rate", value = 0.5, min = 0, max = 2, step = 0.1,
                       ticks = TRUE, width = "100%")
        ),
        div(class = "slider-container",
            sliderInput("eff_condom", "Condom Effectiveness", value = 0.9, min = 0, max = 1, step = 0.05,
                       ticks = TRUE, width = "100%")
        ),
        
        # HIV natural history
        h5("HIV Natural History"),
        sliderInput("dur_acute", "Acute Infection Duration (months)", value = 3, min = 1, max = 12, step = 1,
                   ticks = TRUE, width = "100%"),
        sliderInput("dur_latent", "Latent Infection Duration (years)", value = 10, min = 1, max = 20, step = 1,
                   ticks = TRUE, width = "100%"),
        sliderInput("dur_falling", "Late-stage Duration (years)", value = 3, min = 1, max = 10, step = 1,
                   ticks = TRUE, width = "100%"),
        sliderInput("rel_trans_acute", "Acute Transmissibility Multiplier", value = 6, min = 1, max = 20, step = 1,
                   ticks = TRUE, width = "100%"),
        sliderInput("rel_trans_falling", "Late-stage Transmissibility Multiplier", value = 8, min = 1, max = 20, step = 1,
                   ticks = TRUE, width = "100%"),
        
        # Interventions
        h5("Interventions"),
        checkboxInput("include_testing", "Include HIV Testing", value = TRUE),
        checkboxInput("include_art", "Include ART Treatment", value = TRUE),
        checkboxInput("include_vmmc", "Include VMMC", value = TRUE),
        checkboxInput("include_prep", "Include PrEP", value = FALSE),
        
        # Network parameters
        h5("Sexual Network"),
        selectInput("network_type", "Network Type", 
                   choices = list("Structured Sexual" = "structured",
                                "MSM" = "msm",
                                "Mixed" = "mixed"),
                   selected = "structured"),
        sliderInput("debut_age", "Sexual Debut Age", value = 20, min = 15, max = 30, step = 1,
                   ticks = TRUE, width = "100%"),
        
        # Action buttons
        br(),
        actionButton("run_simulation", "Run Simulation", class = "btn-primary", width = "100%"),
        br(), br(),
        actionButton("reset_params", "Reset Parameters", class = "btn-secondary", width = "100%")
    )
  ),
  
  dashboardBody(
    # Custom CSS for better sidebar appearance
    tags$head(
      tags$style(HTML("
        .sidebar {
          width: 350px !important;
        }
        .main-sidebar {
          width: 350px !important;
        }
        .content-wrapper {
          margin-left: 350px !important;
        }
        .parameter_controls {
          padding: 15px;
        }
        .slider-container {
          margin-bottom: 15px;
        }
        .slider-container .irs {
          margin-top: 10px;
        }
        .nav-tabs {
          border-bottom: 2px solid #dee2e6;
        }
        .nav-tabs .nav-link {
          border: 1px solid transparent;
          border-top-left-radius: 0.25rem;
          border-top-right-radius: 0.25rem;
          margin-bottom: -1px;
        }
        .nav-tabs .nav-link.active {
          color: #495057;
          background-color: #fff;
          border-color: #dee2e6 #dee2e6 #fff;
        }
        .nav-tabs .nav-link:hover {
          border-color: #e9ecef #e9ecef #dee2e6;
        }
      "))
    ),
    
    # Main content area with tab panels
    fluidRow(
      column(width = 12,
        tabsetPanel(
          id = "main_tabs",
          type = "tabs",
          
          # Simulation Tab
          tabPanel("Simulation", icon = icon("play"),
            fluidRow(
              box(title = "Simulation Status", status = "info", solidHeader = TRUE, width = 12,
                  textOutput("simulation_status"),
                  textOutput("progress_text")
              )
            ),
            
            fluidRow(
              box(title = "Quick Results", status = "success", solidHeader = TRUE, width = 6,
                  plotlyOutput("quick_plot", height = "400px")
              ),
              box(title = "Population Summary", status = "info", solidHeader = TRUE, width = 6,
                  tableOutput("population_summary")
              )
            )
          ),
          
          # Results Tab
          tabPanel("Results", icon = icon("chart-line"),
            fluidRow(
              box(title = "HIV Prevalence Over Time", status = "primary", solidHeader = TRUE, width = 12,
                  plotlyOutput("prevalence_plot", height = "500px")
              )
            ),
            
            fluidRow(
              box(title = "HIV Incidence", status = "warning", solidHeader = TRUE, width = 6,
                  plotlyOutput("incidence_plot", height = "400px")
              ),
              box(title = "CD4 Count Distribution", status = "info", solidHeader = TRUE, width = 6,
                  plotlyOutput("cd4_plot", height = "400px")
              )
            ),
            
            fluidRow(
              box(title = "Treatment Coverage", status = "success", solidHeader = TRUE, width = 6,
                  plotlyOutput("treatment_plot", height = "400px")
              ),
              box(title = "Age-specific Prevalence", status = "danger", solidHeader = TRUE, width = 6,
                  plotlyOutput("age_prev_plot", height = "400px")
              )
            ),
            
            fluidRow(
              box(title = "Network Analysis", status = "info", solidHeader = TRUE, width = 12,
                  plotlyOutput("network_plot", height = "400px")
              )
            )
          ),
          
          # Parameters Tab
          tabPanel("Parameters", icon = icon("cog"),
            fluidRow(
              box(title = "Current Parameters", status = "info", solidHeader = TRUE, width = 12,
                  DTOutput("parameters_table")
              )
            ),
            
            fluidRow(
              box(title = "Parameter Sensitivity", status = "warning", solidHeader = TRUE, width = 12,
                  plotlyOutput("sensitivity_plot", height = "500px")
              )
            )
          ),
          
          # About Tab
          tabPanel("About", icon = icon("info"),
            fluidRow(
              box(title = "About STIsim HIV Model", status = "info", solidHeader = TRUE, width = 12,
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
  )
)
