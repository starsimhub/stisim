# HIV Simulation Shiny App - Auto Install and Start Script (PowerShell)
# This script automatically installs all dependencies and starts the app
# Run from the shiny/ directory

Write-Host "🚀 Starting HIV Simulation Shiny App Setup..." -ForegroundColor Blue
Write-Host "==============================================" -ForegroundColor Blue

# Check if we're in the right directory
if (-not (Test-Path "app.R") -or -not (Test-Path "python")) {
    Write-Host "[ERROR] Please run this script from the shiny directory (where app.R is located)" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# Kill any existing R processes on port 3030
Write-Host "[INFO] Checking for existing processes on port 3030..." -ForegroundColor Yellow
$processes = Get-NetTCPConnection -LocalPort 3030 -ErrorAction SilentlyContinue
if ($processes) {
    Write-Host "[INFO] Found existing process on port 3030, killing it..." -ForegroundColor Yellow
    $processes | ForEach-Object { 
        $pid = (Get-Process -Id $_.OwningProcess -ErrorAction SilentlyContinue).Id
        if ($pid) { 
            Stop-Process -Id $pid -Force -ErrorAction SilentlyContinue 
        }
    }
    Start-Sleep -Seconds 2
}

# Check if R is installed
Write-Host "[INFO] Checking R installation..." -ForegroundColor Yellow
try {
    $rVersion = R --version 2>$null
    if ($LASTEXITCODE -ne 0) { throw "R not found" }
    Write-Host "[SUCCESS] R is installed" -ForegroundColor Green
} catch {
    Write-Host "[ERROR] R is not installed. Please install R first from https://cran.r-project.org/" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# Check if Python is installed
Write-Host "[INFO] Checking Python installation..." -ForegroundColor Yellow
try {
    $pythonVersion = python --version 2>$null
    if ($LASTEXITCODE -ne 0) { throw "Python not found" }
    Write-Host "[SUCCESS] Python is installed" -ForegroundColor Green
} catch {
    Write-Host "[ERROR] Python is not installed. Please install Python first from https://python.org/" -ForegroundColor Red
    Read-Host "Press Enter to exit"
    exit 1
}

# Go to project root for virtual environment setup
Write-Host "[INFO] Setting up Python virtual environment..." -ForegroundColor Yellow
Set-Location ..

# Create virtual environment if it doesn't exist
if (-not (Test-Path "venv")) {
    Write-Host "[INFO] Creating Python virtual environment..." -ForegroundColor Yellow
    python -m venv venv
}

# Activate virtual environment
Write-Host "[INFO] Activating virtual environment..." -ForegroundColor Yellow
& "venv\Scripts\Activate.ps1"

# Upgrade pip
Write-Host "[INFO] Upgrading pip..." -ForegroundColor Yellow
python -m pip install --upgrade pip

# Install Python dependencies
Write-Host "[INFO] Installing Python dependencies..." -ForegroundColor Yellow
if (Test-Path "shiny\python\requirements.txt") {
    pip install -r shiny\python\requirements.txt
} else {
    Write-Host "[WARNING] shiny\python\requirements.txt not found, installing basic dependencies..." -ForegroundColor Yellow
    pip install numpy pandas matplotlib seaborn plotly scipy scikit-learn
}

# Install STIsim package in editable mode
Write-Host "[INFO] Installing STIsim package..." -ForegroundColor Yellow
pip install -e .

# Install R dependencies
Write-Host "[INFO] Installing R dependencies..." -ForegroundColor Yellow
Set-Location shiny
$rScript = @"
required_packages <- c('shiny', 'shinydashboard', 'plotly', 'DT', 'reticulate', 'dplyr', 'ggplot2', 'rstudioapi', 'shinyWidgets', 'htmltools', 'jsonlite')
install_if_missing <- function(pkg) { 
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) { 
    cat('Installing', pkg, '...\n'); 
    install.packages(pkg, repos = 'https://cran.r-project.org/') 
  } else { 
    cat(pkg, 'already installed\n') 
  } 
}
for (pkg in required_packages) { install_if_missing(pkg) }
cat('All R packages installed successfully!\n')
"@

$rScript | Rscript -

# Final check - kill any processes that might have started
Write-Host "[INFO] Final cleanup - ensuring port 3030 is free..." -ForegroundColor Yellow
$processes = Get-NetTCPConnection -LocalPort 3030 -ErrorAction SilentlyContinue
if ($processes) {
    $processes | ForEach-Object { 
        $pid = (Get-Process -Id $_.OwningProcess -ErrorAction SilentlyContinue).Id
        if ($pid) { 
            Stop-Process -Id $pid -Force -ErrorAction SilentlyContinue 
        }
    }
}

# Start the app
Write-Host "[INFO] Starting the HIV Simulation Shiny App..." -ForegroundColor Yellow
Write-Host "[SUCCESS] Setup complete! Starting app on http://localhost:3030" -ForegroundColor Green
Write-Host "[INFO] Press Ctrl+C to stop the app" -ForegroundColor Yellow

Rscript app.R