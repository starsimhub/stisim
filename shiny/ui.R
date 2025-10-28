# STIsim HIV Shiny Web App - UI
# User interface definition

ui <- dashboardPage(
  dashboardHeader(
    title = div(
      icon("chart-line", class = "fa-lg"),
      "STIsim HIV Model",
      style = "display: flex; align-items: center; gap: 10px;"
    ),
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    
    # Parameter controls with collapsible sections
    div(id = "parameter_controls",
        h4("Simulation Parameters"),
        
        # Basic simulation parameters
        div(class = "parameter-section",
            h5(class = "parameter-section-header", 
               icon("cog", class = "section-icon"),
               "Basic Parameters",
               icon("chevron-down", class = "collapse-icon")
            ),
            div(class = "parameter-section-content",
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
                )
            )
        ),
        
        # HIV disease parameters
        div(class = "parameter-section",
            h5(class = "parameter-section-header", 
               icon("heartbeat", class = "section-icon"),
               "HIV Parameters",
               icon("chevron-down", class = "collapse-icon")
            ),
            div(class = "parameter-section-content",
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
                )
            )
        ),
        
        # HIV natural history
        div(class = "parameter-section",
            h5(class = "parameter-section-header", 
               icon("clock", class = "section-icon"),
               "HIV Natural History",
               icon("chevron-down", class = "collapse-icon")
            ),
            div(class = "parameter-section-content",
                div(class = "slider-container",
                    sliderInput("dur_acute", "Acute Infection Duration (months)", value = 3, min = 1, max = 12, step = 1,
                               ticks = TRUE, width = "100%")
                ),
                div(class = "slider-container",
                    sliderInput("dur_latent", "Latent Infection Duration (years)", value = 10, min = 1, max = 20, step = 1,
                               ticks = TRUE, width = "100%")
                ),
                div(class = "slider-container",
                    sliderInput("dur_falling", "Late-stage Duration (years)", value = 3, min = 1, max = 10, step = 1,
                               ticks = TRUE, width = "100%")
                ),
                div(class = "slider-container",
                    sliderInput("rel_trans_acute", "Acute Transmissibility Multiplier", value = 6, min = 1, max = 20, step = 1,
                               ticks = TRUE, width = "100%")
                ),
                div(class = "slider-container",
                    sliderInput("rel_trans_falling", "Late-stage Transmissibility Multiplier", value = 8, min = 1, max = 20, step = 1,
                               ticks = TRUE, width = "100%")
                )
            )
        ),
        
        # Interventions
        div(class = "parameter-section",
            h5(class = "parameter-section-header", 
               icon("shield-alt", class = "section-icon"),
               "Interventions",
               icon("chevron-down", class = "collapse-icon")
            ),
            div(class = "parameter-section-content",
                div(class = "checkbox-container",
                    checkboxInput("include_testing", "Include HIV Testing", value = TRUE)
                ),
                div(class = "checkbox-container",
                    checkboxInput("include_art", "Include ART Treatment", value = TRUE)
                ),
                div(class = "checkbox-container",
                    checkboxInput("include_vmmc", "Include VMMC", value = TRUE)
                ),
                div(class = "checkbox-container",
                    checkboxInput("include_prep", "Include PrEP", value = FALSE)
                )
            )
        ),
        
        # Network parameters
        div(class = "parameter-section",
            h5(class = "parameter-section-header", 
               icon("project-diagram", class = "section-icon"),
               "Sexual Network",
               icon("chevron-down", class = "collapse-icon")
            ),
            div(class = "parameter-section-content",
                div(class = "select-container",
                    selectInput("network_type", "Network Type", 
                               choices = list("Structured Sexual" = "structured",
                                            "MSM" = "msm",
                                            "Mixed" = "mixed"),
                               selected = "structured")
                ),
                div(class = "slider-container",
                    sliderInput("debut_age", "Sexual Debut Age", value = 20, min = 15, max = 30, step = 1,
                               ticks = TRUE, width = "100%")
                )
            )
        ),
        
        # Action buttons
        div(class = "action-buttons",
            actionButton("run_simulation", "Run Simulation", class = "btn-primary", width = "90%"),
            br(), br(),
            actionButton("reset_params", "Reset Parameters", class = "btn-secondary", width = "90%")
        )
    )
  ),
  
  dashboardBody(
    # Enhanced CSS for modern, professional appearance
    tags$head(
      tags$style(HTML("
        /* Import Google Fonts */
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
        
        /* Global Styles */
        body {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          min-height: 100vh;
        }
        
        /* Sidebar Styling */
        .sidebar {
          width: 350px !important;
          background: linear-gradient(180deg, #2c3e50 0%, #34495e 100%) !important;
          box-shadow: 2px 0 10px rgba(0,0,0,0.1);
        }
        
        .main-sidebar {
          width: 350px !important;
        }
        
        .content-wrapper {
          margin-left: 350px !important;
          background: transparent;
        }
        
        /* Sidebar Menu Styling */
        .sidebar-menu .treeview-menu > li > a {
          color: #bdc3c7;
          padding: 8px 20px;
          border-radius: 6px;
          margin: 2px 10px;
          transition: all 0.3s ease;
        }
        
        .sidebar-menu .treeview-menu > li > a:hover {
          background: rgba(255,255,255,0.1);
          color: #ffffff;
          transform: translateX(5px);
        }
        
        /* Parameter Controls */
        .parameter_controls {
          padding: 20px;
          background: rgba(255,255,255,0.05);
          border-radius: 10px;
          margin: 10px;
          backdrop-filter: blur(10px);
        }
        
        .parameter_controls h4 {
          color: #ffffff;
          font-weight: 600;
          margin-bottom: 20px;
          text-align: center;
          font-size: 1.2em;
        }
        
        /* Collapsible Parameter Sections */
        .parameter-section {
          margin-bottom: 15px;
          border-radius: 8px;
          background: rgba(255,255,255,0.03);
          border: 1px solid rgba(255,255,255,0.1);
          overflow: hidden;
        }
        
        .parameter-section-header {
          color: #ecf0f1 !important;
          font-weight: 600 !important;
          margin: 0 !important;
          padding: 12px 15px !important;
          background: rgba(255,255,255,0.1) !important;
          cursor: pointer !important;
          display: flex !important;
          align-items: center !important;
          justify-content: space-between !important;
          transition: all 0.3s ease !important;
          font-size: 0.95em !important;
        }
        
        .parameter-section-header:hover {
          background: rgba(255,255,255,0.15) !important;
          color: #ffffff !important;
        }
        
        .parameter-section-header.collapsed {
          background: rgba(255,255,255,0.05) !important;
        }
        
        .section-icon {
          margin-right: 8px;
          color: #3498db;
          font-size: 0.9em;
        }
        
        .collapse-icon {
          transition: transform 0.3s ease;
          color: #bdc3c7;
          font-size: 0.8em;
        }
        
        .parameter-section-header.collapsed .collapse-icon {
          transform: rotate(-90deg);
        }
        
        .parameter-section-content {
          padding: 15px;
          background: rgba(255,255,255,0.02);
          transition: all 0.3s ease;
          max-height: 1000px;
          overflow: hidden;
        }
        
        .parameter-section-content.collapsed {
          max-height: 0;
          padding: 0 15px;
        }
        
        .checkbox-container {
          margin: 8px 0;
          padding: 8px;
          background: rgba(255,255,255,0.05);
          border-radius: 6px;
          border: 1px solid rgba(255,255,255,0.1);
        }
        
        .select-container {
          margin: 8px 0;
          padding: 8px;
          background: rgba(255,255,255,0.05);
          border-radius: 6px;
          border: 1px solid rgba(255,255,255,0.1);
        }
        
        .action-buttons {
          margin-top: 20px;
          padding-top: 15px;
          border-top: 1px solid rgba(255,255,255,0.1);
          text-align: center;
        }
        
        .action-buttons .btn {
          margin: 5px 0;
          min-width: 200px;
          width: 70% !important;
        }
        
        /* Override Shiny's inline styles for action buttons */
        .action-buttons .action-button {
          width: 70% !important;
        }
        
        .action-buttons #run_simulation {
          width: 70% !important;
        }
        
        .action-buttons #reset_params {
          width: 70% !important;
        }
        
        /* Slider Styling */
        .slider-container {
          margin-bottom: 20px;
          background: rgba(255,255,255,0.1);
          padding: 15px;
          border-radius: 8px;
          border: 1px solid rgba(255,255,255,0.2);
        }
        
        .slider-container .irs {
          margin-top: 10px;
        }
        
        .slider-container label {
          color: #ecf0f1 !important;
          font-weight: 500;
          font-size: 0.9em;
        }
        
        /* Custom Slider Colors */
        .irs-bar {
          background: linear-gradient(90deg, #3498db, #2980b9) !important;
        }
        
        .irs-handle {
          background: #ffffff !important;
          border: 3px solid #3498db !important;
          box-shadow: 0 2px 6px rgba(0,0,0,0.3);
        }
        
        .irs-from, .irs-to, .irs-single {
          background: #3498db !important;
          color: white !important;
        }
        
        /* Checkbox Styling */
        .checkbox {
          margin: 10px 0;
        }
        
        .checkbox label {
          color: #ecf0f1 !important;
          font-weight: 400;
        }
        
        /* Select Input Styling */
        .selectize-input {
          background: rgba(255,255,255,0.1) !important;
          border: 1px solid rgba(255,255,255,0.3) !important;
          color: #ecf0f1 !important;
          border-radius: 6px;
        }
        
        .selectize-dropdown {
          background: #34495e !important;
          border: 1px solid rgba(255,255,255,0.2) !important;
        }
        
        /* Button Styling */
        .btn-primary {
          background: linear-gradient(45deg, #3498db, #2980b9) !important;
          border: none !important;
          border-radius: 8px !important;
          font-weight: 600 !important;
          padding: 12px 20px !important;
          box-shadow: 0 4px 15px rgba(52, 152, 219, 0.3) !important;
          transition: all 0.3s ease !important;
        }
        
        .btn-primary:hover {
          transform: translateY(-2px) !important;
          box-shadow: 0 6px 20px rgba(52, 152, 219, 0.4) !important;
        }
        
        .btn-secondary {
          background: linear-gradient(45deg, #95a5a6, #7f8c8d) !important;
          border: none !important;
          border-radius: 8px !important;
          font-weight: 600 !important;
          padding: 12px 20px !important;
          box-shadow: 0 4px 15px rgba(149, 165, 166, 0.3) !important;
          transition: all 0.3s ease !important;
        }
        
        .btn-secondary:hover {
          transform: translateY(-2px) !important;
          box-shadow: 0 6px 20px rgba(149, 165, 166, 0.4) !important;
        }
        
        /* Tab Styling */
        .nav-tabs {
          border-bottom: 3px solid rgba(255,255,255,0.2);
          background: rgba(255,255,255,0.1);
          border-radius: 10px 10px 0 0;
          padding: 5px 5px 0 5px;
          backdrop-filter: blur(10px);
        }
        
        .nav-tabs .nav-link {
          border: none !important;
          border-radius: 8px 8px 0 0 !important;
          margin: 0 2px !important;
          padding: 12px 20px !important;
          color: #ecf0f1 !important;
          font-weight: 500 !important;
          transition: all 0.3s ease !important;
          background: transparent !important;
        }
        
        .nav-tabs .nav-link:hover {
          background: rgba(255,255,255,0.1) !important;
          color: #ffffff !important;
          transform: translateY(-2px);
        }
        
        .nav-tabs .nav-link.active {
          color: #2c3e50 !important;
          background: linear-gradient(135deg, #ffffff, #f8f9fa) !important;
          border: none !important;
          box-shadow: 0 -2px 10px rgba(0,0,0,0.1) !important;
        }
        
        /* Box Styling */
        .box {
          border-radius: 12px !important;
          box-shadow: 0 8px 32px rgba(0,0,0,0.1) !important;
          border: none !important;
          background: rgba(255,255,255,0.95) !important;
          backdrop-filter: blur(10px) !important;
          margin-bottom: 20px !important;
        }
        
        .box-header {
          border-radius: 12px 12px 0 0 !important;
          background: linear-gradient(135deg, #667eea, #764ba2) !important;
          color: white !important;
          padding: 15px 20px !important;
        }
        
        .box-title {
          font-weight: 600 !important;
          font-size: 1.1em !important;
        }
        
        .box-body {
          padding: 20px !important;
        }
        
        /* Status Boxes */
        .box.box-info {
          border-left: 4px solid #3498db !important;
        }
        
        .box.box-success {
          border-left: 4px solid #27ae60 !important;
        }
        
        .box.box-warning {
          border-left: 4px solid #f39c12 !important;
        }
        
        .box.box-danger {
          border-left: 4px solid #e74c3c !important;
        }
        
        .box.box-primary {
          border-left: 4px solid #9b59b6 !important;
        }
        
        /* Table Styling */
        .table {
          border-radius: 8px !important;
          overflow: hidden !important;
        }
        
        .table thead th {
          background: linear-gradient(135deg, #667eea, #764ba2) !important;
          color: white !important;
          font-weight: 600 !important;
          border: none !important;
        }
        
        .table tbody tr:nth-child(even) {
          background: rgba(52, 152, 219, 0.05) !important;
        }
        
        .table tbody tr:hover {
          background: rgba(52, 152, 219, 0.1) !important;
        }
        
        /* Progress Bar Styling */
        .progress {
          height: 8px !important;
          border-radius: 4px !important;
          background: rgba(255,255,255,0.2) !important;
        }
        
        .progress-bar {
          background: linear-gradient(90deg, #3498db, #2980b9) !important;
          border-radius: 4px !important;
        }
        
        /* Well Panel Styling */
        .well {
          background: rgba(255,255,255,0.1) !important;
          border: 1px solid rgba(255,255,255,0.2) !important;
          border-radius: 10px !important;
          backdrop-filter: blur(10px) !important;
        }
        
        /* Header Styling */
        .main-header .navbar {
          background: linear-gradient(135deg, #2c3e50, #34495e) !important;
          box-shadow: 0 2px 10px rgba(0,0,0,0.1) !important;
        }
        
        .main-header .navbar .navbar-brand {
          color: #ffffff !important;
          font-weight: 700 !important;
          font-size: 1.3em !important;
        }
        
        /* Responsive Design */
        @media (max-width: 768px) {
          .sidebar {
            width: 100% !important;
          }
          .content-wrapper {
            margin-left: 0 !important;
          }
        }
        
        /* Animation for smooth transitions */
        * {
          transition: all 0.3s ease;
        }
        
        /* Custom scrollbar */
        ::-webkit-scrollbar {
          width: 8px;
        }
        
        ::-webkit-scrollbar-track {
          background: rgba(255,255,255,0.1);
        }
        
        ::-webkit-scrollbar-thumb {
          background: rgba(255,255,255,0.3);
          border-radius: 4px;
        }
        
        ::-webkit-scrollbar-thumb:hover {
          background: rgba(255,255,255,0.5);
        }
      ")),
      
      # JavaScript for collapsible sections and button width fix
      tags$script(HTML("
        $(document).ready(function() {
          // Add click handlers to all parameter section headers
          $('.parameter-section-header').click(function() {
            var content = $(this).next('.parameter-section-content');
            var icon = $(this).find('.collapse-icon');
            
            // Toggle the collapsed class
            $(this).toggleClass('collapsed');
            content.toggleClass('collapsed');
            
            // Rotate the chevron icon
            if ($(this).hasClass('collapsed')) {
              icon.css('transform', 'rotate(-90deg)');
            } else {
              icon.css('transform', 'rotate(0deg)');
            }
          });
          
          // Initialize all sections as collapsed
          $('.parameter-section-header').each(function() {
            var content = $(this).next('.parameter-section-content');
            var icon = $(this).find('.collapse-icon');
            
            // Set initial state as collapsed
            $(this).addClass('collapsed');
            content.addClass('collapsed');
            icon.css('transform', 'rotate(-90deg)');
          });
          
          // Fix button widths by overriding inline styles
          function fixButtonWidths() {
            $('#run_simulation, #reset_params').css('width', '70%');
          }
          
          // Apply immediately and after any Shiny updates
          fixButtonWidths();
          $(document).on('shiny:value', fixButtonWidths);
          $(document).on('shiny:outputinvalidated', fixButtonWidths);
        });
      "))
    ),
    
    # Main content area with tab panels
    fluidRow(
      column(width = 12,
        tabsetPanel(
          id = "main_tabs",
          type = "tabs",
          
          # Simulation & Results Tab
          tabPanel("Simulation & Results", icon = icon("play-circle"),
            div(style = "padding: 20px;",
              # Simulation Status and Population Summary
              fluidRow(
                box(title = div(icon("info-circle"), "Simulation Status"), status = "info", solidHeader = TRUE, width = 6,
                    div(style = "font-size: 16px; padding: 10px;",
                        textOutput("simulation_status")
                    ),
                    div(style = "margin-top: 10px;",
                        textOutput("progress_text")
                    )
                ),
                box(title = div(icon("users"), "Population Summary"), status = "info", solidHeader = TRUE, width = 6,
                    div(style = "overflow-x: auto;",
                        tableOutput("population_summary")
                    )
                )
              ),
              
              # Main Results Section
              fluidRow(
                box(title = div(icon("heartbeat"), "HIV Prevalence Over Time"), status = "primary", solidHeader = TRUE, width = 12,
                    div(style = "text-align: center;",
                        plotlyOutput("prevalence_plot", height = "500px")
                    )
                )
              ),
              
              # Secondary Results Row 1
              fluidRow(
                box(title = div(icon("line-chart"), "HIV Incidence"), status = "warning", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("incidence_plot", height = "400px")
                    )
                ),
                box(title = div(icon("chart-area"), "CD4 Count Distribution"), status = "info", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("cd4_plot", height = "400px")
                    )
                )
              ),
              
              # Secondary Results Row 2
              fluidRow(
                box(title = div(icon("pills"), "Treatment Coverage"), status = "success", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("treatment_plot", height = "400px")
                    )
                ),
                box(title = div(icon("calendar-alt"), "Age-specific Prevalence"), status = "danger", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("age_prev_plot", height = "400px")
                    )
                )
              ),
              
              # Network Analysis
              fluidRow(
                box(title = div(icon("project-diagram"), "Network Analysis"), status = "info", solidHeader = TRUE, width = 12,
                    div(style = "text-align: center;",
                        plotlyOutput("network_analysis_plot", height = "400px")
                    )
                )
              )
            )
          ),
          
          # Parameters Tab
          tabPanel("Parameters", icon = icon("cogs"),
            div(style = "padding: 20px;",
              fluidRow(
                box(title = div(icon("list-alt"), "Current Parameters"), status = "info", solidHeader = TRUE, width = 12,
                    div(style = "overflow-x: auto;",
                        DTOutput("parameters_table")
                    )
                )
              ),
              
              fluidRow(
                box(title = div(icon("chart-pie"), "Parameter Sensitivity Analysis"), status = "warning", solidHeader = TRUE, width = 12,
                    div(style = "text-align: center;",
                        plotlyOutput("sensitivity_plot", height = "500px")
                    )
                )
              )
            )
          ),
          
          # Advanced Analytics Tab
          tabPanel("Advanced Analytics", icon = icon("chart-line"),
            div(style = "padding: 20px;",
              # Time Series Analysis
              fluidRow(
                box(title = div(icon("wave-square"), "HIV Prevalence Trends by Age Group"), status = "primary", solidHeader = TRUE, width = 12,
                    div(style = "text-align: center;",
                        plotlyOutput("age_prevalence_plot", height = "500px")
                    )
                )
              ),
              
              # Heatmap and Distribution Analysis
              fluidRow(
                box(title = div(icon("fire"), "Transmission Risk Heatmap"), status = "danger", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("transmission_heatmap", height = "400px")
                    )
                ),
                box(title = div(icon("chart-bar"), "CD4 Count Distribution"), status = "success", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("cd4_distribution", height = "400px")
                    )
                )
              ),
              
              # Network and Intervention Analysis
              fluidRow(
                box(title = div(icon("project-diagram"), "Sexual Network Visualization"), status = "warning", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("network_visualization_plot", height = "400px")
                    )
                ),
                box(title = div(icon("shield-alt"), "Intervention Impact Analysis"), status = "info", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("intervention_impact", height = "400px")
                    )
                )
              ),
              
              # Advanced Statistical Plots
              fluidRow(
                box(title = div(icon("chart-area"), "Cumulative Infections Over Time"), status = "primary", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("cumulative_infections", height = "400px")
                    )
                ),
                box(title = div(icon("bar-chart"), "Prevalence vs Incidence Correlation"), status = "success", solidHeader = TRUE, width = 6,
                    div(style = "text-align: center;",
                        plotlyOutput("prevalence_incidence_scatter", height = "400px")
                    )
                )
              ),
              
              # Survival and Treatment Analysis
              fluidRow(
                box(title = div(icon("heartbeat"), "Survival Analysis by Treatment Status"), status = "danger", solidHeader = TRUE, width = 12,
                    div(style = "text-align: center;",
                        plotlyOutput("survival_analysis", height = "500px")
                    )
                )
              )
            )
          ),
          
          # About Tab
          tabPanel("About", icon = icon("info-circle"),
            div(style = "padding: 20px;",
              fluidRow(
                box(title = div(icon("microscope"), "About STIsim HIV Model"), status = "info", solidHeader = TRUE, width = 12,
                    div(style = "font-size: 16px; line-height: 1.6;",
                        h4(style = "color: #2c3e50; margin-bottom: 20px;", 
                           icon("chart-line", style = "margin-right: 10px;"), 
                           "STIsim HIV Web Application"),
                        
                        p(style = "font-size: 18px; color: #34495e; margin-bottom: 20px;",
                          "This web application provides an interactive interface for the STIsim HIV modeling framework."),
                        
                        p(style = "font-size: 16px; color: #7f8c8d; margin-bottom: 30px;",
                          "STIsim is an agent-based modeling framework for simulating sexually transmitted diseases, built on the Starsim architecture."),
                        
                        fluidRow(
                          column(6,
                            h5(style = "color: #2c3e50; margin-bottom: 15px; font-weight: 600;",
                               icon("star", style = "margin-right: 8px;"), "Key Features:"),
                            tags$ul(style = "list-style: none; padding: 0;",
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("users"), " Agent-based HIV simulation"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("sliders-h"), " Configurable disease parameters"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("project-diagram"), " Sexual network modeling"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("shield-alt"), " Intervention modeling"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("chart-bar"), " Interactive visualizations"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(52, 152, 219, 0.1); border-radius: 5px;",
                                     icon("chart-pie"), " Parameter sensitivity analysis")
                            )
                          ),
                          column(6,
                            h5(style = "color: #2c3e50; margin-bottom: 15px; font-weight: 600;",
                               icon("cogs", style = "margin-right: 8px;"), "Model Components:"),
                            tags$ul(style = "list-style: none; padding: 0;",
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(46, 204, 113, 0.1); border-radius: 5px;",
                                     icon("heartbeat"), " HIV natural history (acute, latent, late-stage)"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(46, 204, 113, 0.1); border-radius: 5px;",
                                     icon("network-wired"), " Sexual transmission networks"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(46, 204, 113, 0.1); border-radius: 5px;",
                                     icon("user-friends"), " Demographic processes"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(46, 204, 113, 0.1); border-radius: 5px;",
                                     icon("pills"), " HIV testing and treatment"),
                              tags$li(style = "margin: 8px 0; padding: 8px; background: rgba(46, 204, 113, 0.1); border-radius: 5px;",
                                     icon("shield-virus"), " Prevention interventions")
                            )
                          )
                        ),
                        
                        hr(style = "border: 1px solid #ecf0f1; margin: 30px 0;"),
                        
                        fluidRow(
                          column(6,
                            h5(style = "color: #2c3e50; margin-bottom: 15px; font-weight: 600;",
                               icon("envelope", style = "margin-right: 8px;"), "Contact & Support:"),
                            p(style = "color: #7f8c8d; margin-bottom: 10px;",
                              "For questions or support, please visit:"),
                            tags$ul(style = "list-style: none; padding: 0;",
                              tags$li(style = "margin: 5px 0;",
                                     icon("github"), " STIsim GitHub repository"),
                              tags$li(style = "margin: 5px 0;",
                                     icon("at"), " info@starsim.org")
                            )
                          ),
                          column(6,
                            h5(style = "color: #2c3e50; margin-bottom: 15px; font-weight: 600;",
                               icon("code-branch", style = "margin-right: 8px;"), "Version Information:"),
                            div(style = "background: rgba(52, 152, 219, 0.1); padding: 15px; border-radius: 8px;",
                              p(style = "margin: 5px 0; font-weight: 600; color: #2c3e50;",
                                "STIsim v1.4"),
                              p(style = "margin: 5px 0; color: #7f8c8d;",
                                "Shiny Web Interface v1.0"),
                              p(style = "margin: 5px 0; color: #7f8c8d; font-size: 14px;",
                                "© 2024-2025 by IDM")
                            )
                          )
                        )
                    )
                )
              )
            )
          )
        )
      )
    )
  )
)
