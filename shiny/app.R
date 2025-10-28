# STIsim HIV Shiny Web App
# Main application file

# Load required libraries
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(reticulate)
library(dplyr)
library(ggplot2)

# Set up Python environment with error handling
python_available <- FALSE
tryCatch({
  # Use the virtual environment Python where stisim is installed
  use_python("../venv/bin/python", required = FALSE)
  source_python("python/hiv_simulator.py")
  python_available <<- TRUE
  cat("STIsim Python environment configured successfully\n")
}, error = function(e) {
  cat("Python integration not available:", e$message, "\n")
  cat("Running in simulation-only mode\n")
})

# Source UI and server components
source("ui.R")
source("server.R")

# Run the application on port 3838
shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
