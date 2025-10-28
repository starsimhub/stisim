# STIsim HIV Shiny Web App - Dependencies

This document lists all the dependencies required for the STIsim HIV Shiny Web App.

## R Dependencies

The following R packages are required and can be installed using `install_all_dependencies.R`:

### Core Shiny Packages
- **shiny** - Web application framework
- **shinydashboard** - Dashboard layout for Shiny
- **plotly** - Interactive plotting
- **DT** - Data tables
- **reticulate** - Python integration

### Data Manipulation & Visualization
- **dplyr** - Data manipulation
- **ggplot2** - Static plotting
- **htmltools** - HTML generation
- **jsonlite** - JSON handling

### Additional UI Components
- **shinyWidgets** - Enhanced UI widgets
- **rstudioapi** - RStudio integration

## Python Dependencies

The following Python packages are required and can be installed using `python/requirements.txt`:

### Core STIsim Dependencies
- **starsim>=3.0.1** - Agent-based simulation framework
- **sciris>=3.1.6** - Scientific computing utilities
- **pandas>=2.0.0** - Data manipulation
- **scipy** - Scientific computing
- **numba** - JIT compilation
- **networkx** - Network analysis

### STIsim Specific
- **optuna** - Hyperparameter optimization
- **requests** - HTTP library

### Web App Visualization
- **plotly>=5.0.0** - Interactive plotting
- **numpy>=1.20.0** - Numerical computing
- **matplotlib>=3.5.0** - Static plotting

### Additional Scientific Computing
- **seaborn** - Statistical visualization
- **openpyxl** - Excel file handling
- **xlsxwriter** - Excel writing
- **psutil** - System monitoring
- **dill** - Serialization
- **zstandard** - Compression
- **multiprocess** - Parallel processing
- **jsonpickle** - JSON serialization
- **setuptools** - Package building
- **gitpython** - Git integration
- **jellyfish** - String matching
- **memory_profiler** - Memory profiling
- **line_profiler** - Line profiling

## Installation

### Automatic Installation
Run the comprehensive installer:
```bash
cd shiny
Rscript install_all_dependencies.R
```

### Manual Installation

#### R Packages
```r
# Install R packages
source("install_dependencies.R")
```

#### Python Packages
```bash
# Install Python packages
pip install -r python/requirements.txt

# Install STIsim in editable mode
pip install -e ..
```

## Verification

After installation, verify everything is working:

1. **R Dependencies**: All packages should load without errors
2. **Python Integration**: STIsim should be importable
3. **App Launch**: `Rscript app.R` should start the server successfully

## System Requirements

- **R**: Version 4.0 or higher
- **Python**: Version 3.9 or higher
- **Operating System**: Windows, macOS, or Linux
- **Memory**: At least 4GB RAM (8GB+ recommended for large simulations)
- **Storage**: At least 2GB free space

## Troubleshooting

### Common Issues

1. **Python Integration Fails**
   - Ensure Python virtual environment is activated
   - Check that reticulate is using the correct Python path

2. **STIsim Import Error**
   - Verify STIsim is installed in editable mode: `pip install -e ..`
   - Check that all dependencies are installed

3. **Memory Issues**
   - Reduce population size in simulation parameters
   - Close other applications to free up memory

4. **Port Already in Use**
   - Kill existing processes: `lsof -ti:3030 | xargs kill -9`
   - Or change port in `app.R`

### Getting Help

- Check the terminal output for specific error messages
- Verify all dependencies are installed correctly
- Ensure Python virtual environment is properly configured
