#!/bin/bash

# HIV Simulation Shiny App - Auto Install and Start Script
# This script automatically installs all dependencies and starts the app
# Run from the shiny/ directory

set -e  # Exit on any error

echo "🚀 Starting HIV Simulation Shiny App Setup..."
echo "=============================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if we're in the right directory
if [ ! -f "app.R" ] || [ ! -d "python" ]; then
    print_error "Please run this script from the shiny/ directory (where app.R is located)"
    exit 1
fi

# Kill any existing R processes on port 3030
print_status "Checking for existing processes on port 3030..."
if lsof -i :3030 >/dev/null 2>&1; then
    print_warning "Found existing process on port 3030, killing it..."
    lsof -ti :3030 | xargs kill -9 2>/dev/null || true
    sleep 2
fi

# Check if R is installed
print_status "Checking R installation..."
if ! command -v R &> /dev/null; then
    print_error "R is not installed. Please install R first:"
    print_error "  - macOS: brew install r"
    print_error "  - Ubuntu: sudo apt-get install r-base"
    print_error "  - Windows: Download from https://cran.r-project.org/"
    exit 1
fi

# Check if Python is installed
print_status "Checking Python installation..."
if ! command -v python3 &> /dev/null; then
    print_error "Python 3 is not installed. Please install Python 3 first."
    exit 1
fi

# Go to project root for virtual environment setup
print_status "Setting up Python virtual environment..."
cd ..

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    print_status "Creating Python virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
print_status "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
print_status "Upgrading pip..."
pip install --upgrade pip

# Install Python dependencies
print_status "Installing Python dependencies..."
if [ -f "shiny/python/requirements.txt" ]; then
    pip install -r shiny/python/requirements.txt
else
    print_warning "shiny/python/requirements.txt not found, installing basic dependencies..."
    pip install numpy pandas matplotlib seaborn plotly scipy scikit-learn
fi

# Install STIsim package in editable mode
print_status "Installing STIsim package..."
pip install -e .

# Install R dependencies
print_status "Installing R dependencies..."
cd shiny
Rscript -e "
# Install required R packages
required_packages <- c(
  'shiny',
  'shinydashboard', 
  'plotly',
  'DT',
  'reticulate',
  'dplyr',
  'ggplot2',
  'rstudioapi',
  'shinyWidgets',
  'htmltools',
  'jsonlite'
)

# Function to install packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat('Installing', pkg, '...\n')
    install.packages(pkg, repos = 'https://cran.r-project.org/')
  } else {
    cat(pkg, 'already installed\n')
  }
}

# Install all packages
for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat('All R packages installed successfully!\n')
"

# Final check - kill any processes that might have started
print_status "Final cleanup - ensuring port 3030 is free..."
if lsof -i :3030 >/dev/null 2>&1; then
    lsof -ti :3030 | xargs kill -9 2>/dev/null || true
    sleep 2
fi

# Start the app
print_status "Starting the HIV Simulation Shiny App..."
print_success "Setup complete! Starting app on http://localhost:3030"
print_status "Press Ctrl+C to stop the app"

Rscript app.R