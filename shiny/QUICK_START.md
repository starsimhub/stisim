# 🚀 Quick Start - Shiny App

## One-Command Setup and Launch

**Run the startup scripts from the `shiny/` directory:**

### For macOS/Linux:
```bash
cd shiny
./start_app.sh
```

### For Windows:
```cmd
cd shiny
start_app.bat
```

### For Windows PowerShell:
```powershell
cd shiny
.\start_app.ps1
```

## What the Scripts Do

The startup scripts automatically:

1. ✅ **Check Prerequisites** - Verify R and Python are installed
2. ✅ **Kill Existing Processes** - Clean up any running instances on port 3030
3. ✅ **Setup Python Environment** - Create/activate virtual environment in project root
4. ✅ **Install Dependencies** - Install all R and Python packages
5. ✅ **Install STIsim** - Install the simulation package in editable mode
6. ✅ **Start the App** - Launch the Shiny app on http://localhost:3030

## Alternative: Run from Project Root

You can also run the scripts from the project root directory:

```bash
# From project root
./start_app.sh
```

## Manual Steps (if needed)

If the automated scripts don't work, you can run these commands manually:

### 1. Install R Dependencies
```r
# Run in R console or Rscript
source("install_dependencies.R")
```

### 2. Install Python Dependencies
```bash
# From project root
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install packages
pip install -r shiny/python/requirements.txt
pip install -e .  # Install STIsim in editable mode
```

### 3. Start the App
```bash
# From shiny/ directory
Rscript app.R
```

## Troubleshooting

### Port Already in Use
```bash
# Kill processes on port 3030
lsof -ti :3030 | xargs kill -9  # macOS/Linux
# OR
netstat -ano | findstr :3030   # Windows, then taskkill /PID <PID> /F
```

### Missing Dependencies
- **R**: Install from https://cran.r-project.org/
- **Python**: Install from https://python.org/
- **R packages**: The script installs them automatically
- **Python packages**: The script installs them automatically

## App Features

Once running, you'll have access to:
- 🎛️ **Interactive Parameter Controls** - Sliders for all simulation parameters
- 📊 **Real-time Visualizations** - Interactive plots with Plotly
- 🔬 **HIV Simulation** - Powered by STIsim Python package
- 📈 **Multiple Analysis Views** - Prevalence, incidence, and detailed results
- 🎨 **Modern UI** - Clean, responsive dashboard interface

## Access the App

Open your web browser and go to: **http://localhost:3030**
