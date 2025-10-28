# STIsim HIV Shiny Web App

This is a Shiny web application that provides an interactive interface for the STIsim HIV modeling framework.

## Features

- **Interactive HIV Simulation**: Run agent-based HIV simulations with customizable parameters
- **Real-time Visualization**: View HIV prevalence, incidence, and other key metrics
- **Parameter Control**: Adjust disease parameters, intervention settings, and network configurations
- **Multiple Plot Types**: Prevalence trends, age-specific analysis, CD4 distributions, and more
- **Network Analysis**: Visualize sexual network characteristics and partner distributions

## Installation

### Prerequisites

1. **R and RStudio** (or R command line)
2. **Python 3.9+** with pip
3. **STIsim** installed in Python environment

### Automatic Installation (Recommended)

Run the comprehensive installer:

```bash
cd shiny
Rscript install_all_dependencies.R
```

This will install all R packages, Python dependencies, and STIsim automatically.

### Manual Installation

#### R Dependencies

Install required R packages:

```r
source("install_dependencies.R")
```

Or manually:
```r
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "reticulate", "dplyr", "ggplot2", "shinyWidgets", "htmltools", "jsonlite", "rstudioapi"))
```

#### Python Dependencies

Install Python dependencies:

```bash
pip install -r python/requirements.txt
```

#### STIsim Installation

Install STIsim in editable mode:

```bash
pip install -e ..
```

## Usage

### Running the App

1. **From Command Line** (Recommended):
   ```bash
   cd shiny
   Rscript app.R
   ```
   The app will be available at http://localhost:3030

2. **From RStudio**: Open `app.R` and click "Run App"

3. **From R Console**:
   ```r
   setwd("path/to/shiny")
   source("app.R")
   ```

4. **Using the comprehensive installer**:
   ```bash
   cd shiny
   Rscript install_all_dependencies.R
   ```

### Using the Interface

1. **Set Parameters**: Use the sidebar controls to configure simulation parameters
2. **Run Simulation**: Click "Run Simulation" to start the model
3. **View Results**: Navigate through the Results tab to see various visualizations
4. **Analyze Parameters**: Check the Parameters tab for detailed parameter information

## Parameters

### Simulation Parameters
- **Population Size**: Number of agents in the simulation (100-50,000)
- **Start Year**: Beginning of simulation period
- **Duration**: Length of simulation in years

### HIV Disease Parameters
- **Initial Prevalence**: Starting HIV prevalence in the population
- **Transmission Rates**: Male-to-female, mother-to-child, and male-to-male rates
- **Natural History**: Duration of acute, latent, and late-stage infection
- **Transmissibility**: Multipliers for different disease stages

### Interventions
- **HIV Testing**: Voluntary testing with configurable coverage
- **ART Treatment**: Antiretroviral therapy with time-varying coverage
- **VMMC**: Voluntary medical male circumcision
- **PrEP**: Pre-exposure prophylaxis

### Sexual Networks
- **Network Type**: Structured sexual, MSM, or mixed networks
- **Sexual Debut**: Age of sexual debut

## Output Visualizations

### Main Results
- **HIV Prevalence**: Time series of HIV prevalence
- **HIV Incidence**: New infections over time
- **CD4 Distribution**: Distribution of CD4 counts among infected individuals
- **Treatment Coverage**: ART coverage over time

### Advanced Analysis
- **Age-specific Prevalence**: HIV prevalence by age group
- **Network Statistics**: Partner distribution and network characteristics
- **Parameter Sensitivity**: Visualization of parameter values

## File Structure

```
shiny/
├── app.R                          # Main Shiny application
├── ui.R                           # User interface definition
├── server.R                       # Server logic and reactive functions
├── README.md                      # This file
├── ARCHITECTURE.md                # Technical architecture documentation
├── DEPENDENCIES.md                # Dependency documentation
├── install_dependencies.R         # R package installer
├── install_all_dependencies.R     # Comprehensive installer
└── python/
    ├── hiv_simulator.py           # Python backend functions
    └── requirements.txt           # Python dependencies
```

## Troubleshooting

### Common Issues

1. **Python Not Found**: Ensure Python is installed and accessible from R
2. **STIsim Import Error**: Verify STIsim is installed in the correct Python environment
3. **Plot Not Displaying**: Check browser console for JavaScript errors
4. **Simulation Fails**: Verify parameter values are within valid ranges

### Debug Mode

To run in debug mode:

```r
options(shiny.error = browser)
shiny::runApp(display.mode = "showcase")
```

## Model Details

This application uses the STIsim HIV model, which includes:

- **Agent-based simulation** with individual-level modeling
- **Sexual network dynamics** with age-structured partner matching
- **HIV natural history** with acute, latent, and late-stage progression
- **Intervention modeling** for testing, treatment, and prevention
- **Demographic processes** including births, deaths, and aging

## Contact

For questions or support:
- STIsim GitHub: https://github.com/starsimhub/stisim
- Email: info@starsim.org

## License

This application is provided under the MIT License, consistent with the STIsim framework.
