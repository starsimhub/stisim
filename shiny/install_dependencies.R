# STIsim HIV Shiny Web App - Dependency Installer
# Install all required R packages

# List of required packages
required_packages <- c(
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

# Function to install packages if not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# Install all required packages
cat("Installing required R packages...\n")
for (package in required_packages) {
  cat(paste("Installing", package, "...\n"))
  install_if_missing(package)
}

cat("\nAll R packages installed successfully!\n")
cat("Next steps:\n")
cat("1. Install Python dependencies: pip install -r python/requirements.txt\n")
cat("2. Install STIsim: pip install stisim\n")
cat("3. Run the app: Rscript run_app.R\n")
