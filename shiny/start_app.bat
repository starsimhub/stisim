@echo off
REM HIV Simulation Shiny App - Auto Install and Start Script (Windows)
REM This script automatically installs all dependencies and starts the app
REM Run from the shiny/ directory

echo 🚀 Starting HIV Simulation Shiny App Setup...
echo ==============================================

REM Check if we're in the right directory
if not exist "app.R" (
    echo [ERROR] Please run this script from the shiny directory (where app.R is located)
    pause
    exit /b 1
)

if not exist "python" (
    echo [ERROR] Python directory not found
    pause
    exit /b 1
)

REM Kill any existing R processes on port 3030
echo [INFO] Checking for existing processes on port 3030...
for /f "tokens=5" %%a in ('netstat -ano ^| findstr :3030') do (
    echo [INFO] Killing process %%a on port 3030
    taskkill /PID %%a /F >nul 2>&1
)

REM Check if R is installed
echo [INFO] Checking R installation...
where R >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] R is not installed. Please install R first from https://cran.r-project.org/
    pause
    exit /b 1
)

REM Check if Python is installed
echo [INFO] Checking Python installation...
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed. Please install Python first from https://python.org/
    pause
    exit /b 1
)

REM Go to project root for virtual environment setup
echo [INFO] Setting up Python virtual environment...
cd ..

REM Create virtual environment if it doesn't exist
if not exist "venv" (
    echo [INFO] Creating Python virtual environment...
    python -m venv venv
)

REM Activate virtual environment
echo [INFO] Activating virtual environment...
call venv\Scripts\activate.bat

REM Upgrade pip
echo [INFO] Upgrading pip...
python -m pip install --upgrade pip

REM Install Python dependencies
echo [INFO] Installing Python dependencies...
if exist "shiny\python\requirements.txt" (
    pip install -r shiny\python\requirements.txt
) else (
    echo [WARNING] shiny\python\requirements.txt not found, installing basic dependencies...
    pip install numpy pandas matplotlib seaborn plotly scipy scikit-learn
)

REM Install STIsim package in editable mode
echo [INFO] Installing STIsim package...
pip install -e .

REM Install R dependencies
echo [INFO] Installing R dependencies...
cd shiny
Rscript -e "required_packages <- c('shiny', 'shinydashboard', 'plotly', 'DT', 'reticulate', 'dplyr', 'ggplot2', 'rstudioapi', 'shinyWidgets', 'htmltools', 'jsonlite'); install_if_missing <- function(pkg) { if (!require(pkg, character.only = TRUE, quietly = TRUE)) { cat('Installing', pkg, '...\n'); install.packages(pkg, repos = 'https://cran.r-project.org/') } else { cat(pkg, 'already installed\n') } }; for (pkg in required_packages) { install_if_missing(pkg) }; cat('All R packages installed successfully!\n')"

REM Final check - kill any processes that might have started
echo [INFO] Final cleanup - ensuring port 3030 is free...
for /f "tokens=5" %%a in ('netstat -ano ^| findstr :3030') do (
    taskkill /PID %%a /F >nul 2>&1
)

REM Start the app
echo [INFO] Starting the HIV Simulation Shiny App...
echo [SUCCESS] Setup complete! Starting app on http://localhost:3030
echo [INFO] Press Ctrl+C to stop the app

Rscript app.R