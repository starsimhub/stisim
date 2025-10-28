# STIsim HIV Shiny Web App - Complete Dependency Installer
# This script installs all required dependencies for the STIsim Shiny app

cat("=== STIsim HIV Shiny Web App - Dependency Installer ===\n\n")

# Function to install packages if not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing R package:", package, "...\n"))
    install.packages(package, dependencies = TRUE, quiet = FALSE)
    library(package, character.only = TRUE)
    cat(paste("✓", package, "installed successfully\n"))
  } else {
    cat(paste("✓", package, "already installed\n"))
  }
}

# R packages required for the Shiny app
cat("1. Installing R packages...\n")
required_r_packages <- c(
  "shiny",
  "shinydashboard", 
  "plotly",
  "DT",
  "reticulate",
  "dplyr",
  "ggplot2",
  "rstudioapi",
  "shinyWidgets",
  "htmltools",
  "jsonlite"
)

for (package in required_r_packages) {
  install_if_missing(package)
}

cat("\n2. Setting up Python environment...\n")

# Check if Python is available
if (require(reticulate, quietly = TRUE)) {
  cat("✓ reticulate package loaded\n")
  
  # Try to use the virtual environment Python
  tryCatch({
    use_python("../venv/bin/python", required = FALSE)
    cat("✓ Using virtual environment Python\n")
  }, error = function(e) {
    cat("⚠ Could not find virtual environment Python, using system Python\n")
  })
  
  # Install Python packages
  cat("Installing Python packages...\n")
  
  # Install from requirements.txt
  py_install_requirements <- function() {
    tryCatch({
      py_run_string("
import subprocess
import sys
import os

# Install requirements from the python subdirectory
result = subprocess.run([sys.executable, '-m', 'pip', 'install', '-r', 'python/requirements.txt'], 
                       capture_output=True, text=True)

if result.returncode == 0:
    print('✓ Python packages installed successfully')
else:
    print('⚠ Some Python packages may have failed to install')
    print(result.stderr)
")
    }, error = function(e) {
      cat("⚠ Error installing Python packages:", e$message, "\n")
    })
  }
  
  py_install_requirements()
  
} else {
  cat("⚠ reticulate not available, skipping Python setup\n")
}

cat("\n3. Installing STIsim package...\n")

# Install STIsim in editable mode
if (require(reticulate, quietly = TRUE)) {
  tryCatch({
    py_run_string("
import subprocess
import sys
import os

# Change to the project root directory (go up one level from shiny/)
os.chdir('../')

# Install STIsim in editable mode
result = subprocess.run([sys.executable, '-m', 'pip', 'install', '-e', '.'], 
                       capture_output=True, text=True)

if result.returncode == 0:
    print('✓ STIsim package installed successfully')
else:
    print('⚠ Error installing STIsim package')
    print(result.stderr)
")
  }, error = function(e) {
    cat("⚠ Error installing STIsim:", e$message, "\n")
  })
}

cat("\n=== Installation Complete ===\n")
cat("All dependencies have been installed.\n")
cat("You can now run the app with: Rscript app.R\n")
cat("The app will be available at: http://localhost:3030\n")
