# STIsim HIV Shiny Web App - Architecture Documentation

## Overview

The STIsim HIV Shiny Web App is a sophisticated web application that provides an interactive interface for running HIV transmission simulations using the STIsim modeling framework. The application combines R/Shiny for the user interface with Python/STIsim for the simulation engine.

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    STIsim HIV Shiny Web App                    │
├─────────────────────────────────────────────────────────────────┤
│  Frontend (R/Shiny)           │  Backend (Python/STIsim)       │
│  ┌─────────────────────────┐   │  ┌─────────────────────────┐   │
│  │     User Interface      │   │  │   Simulation Engine     │   │
│  │   - Sliders/Controls    │   │  │   - STIsim Framework    │   │
│  │   - Interactive Plots   │   │  │   - HIV Models          │   │
│  │   - Dashboard Layout    │   │  │   - Network Models      │   │
│  └─────────────────────────┘   │  │   - Interventions       │   │
│  ┌─────────────────────────┐   │  └─────────────────────────┘   │
│  │     R Server Logic      │   │  ┌─────────────────────────┐   │
│  │   - Parameter Handling  │   │  │   Data Processing       │   │
│  │   - Plot Generation     │   │  │   - Result Extraction   │   │
│  │   - State Management    │   │  │   - Data Formatting     │   │
│  └─────────────────────────┘   │  └─────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                    ┌─────────────────┐
                    │   Reticulate    │
                    │   Integration   │
                    │   Layer         │
                    └─────────────────┘
```

## Component Breakdown

### 1. Frontend Layer (R/Shiny)

#### **User Interface (`ui.R`)**
- **Dashboard Layout**: Uses `shinydashboard` for professional appearance
- **Parameter Controls**: Interactive sliders for all simulation parameters
- **Visualization Areas**: Multiple tabs for different types of plots and results
- **Responsive Design**: Wide sidebar (350px) for better parameter visibility

**Key UI Components:**
```r
# Main dashboard structure
ui <- dashboardPage(
  dashboardHeader(title = "STIsim HIV Model"),
  dashboardSidebar(width = 350, ...),  # Wide sidebar for parameters
  dashboardBody(...)                   # Main content area
)

# Parameter controls with sliders
sliderInput("n_agents", "Population Size", 
           value = 500, min = 100, max = 50000, step = 100)
```

#### **Server Logic (`server.R`)**
- **Parameter Collection**: Gathers all user inputs from sliders
- **Simulation Orchestration**: Manages the simulation workflow
- **Result Processing**: Handles data from Python and creates plots
- **State Management**: Tracks simulation status and progress

**Key Server Functions:**
```r
# Simulation execution
observeEvent(input$run_simulation, {
  # Collect parameters
  params <- list(
    n_agents = as.integer(input$n_agents),
    start = input$start_year,
    # ... more parameters
  )
  
  # Call Python simulation
  if (python_available) {
    result <- py$run_hiv_simulation(params)
  } else {
    result <- generate_mock_data(params)
  }
})
```

### 2. Integration Layer (Reticulate)

#### **Python Environment Setup (`app.R`)**
The `reticulate` package provides seamless R-Python integration:

```r
# Set up Python environment
use_python("../venv/bin/python", required = FALSE)
source_python("python/hiv_simulator.py")
python_available <<- TRUE
```

**How Reticulate Works:**
1. **`use_python()`**: Specifies which Python interpreter to use
2. **`source_python()`**: Loads Python file and makes functions available via `py$`
3. **`py$function_name()`**: Calls Python functions from R
4. **Automatic Data Conversion**: Converts R objects to Python and vice versa

#### **Data Flow Between R and Python**
```
R Parameters → Python Dictionary → STIsim Simulation → Python Results → R Objects
```

### 3. Backend Layer (Python/STIsim)

#### **Simulation Engine (`python/hiv_simulator.py`)**

**Main Function:**
```python
def run_hiv_simulation(params):
    """Run HIV simulation with given parameters"""
    # Create HIV disease module
    hiv = sti.HIV(
        init_prev=params['init_prev'],
        beta_m2f=params['beta_m2f'],
        # ... more parameters
    )
    
    # Create networks
    networks = []
    if params['network_type'] == 'structured':
        prior_partners = sti.PriorPartners()
        sexual = sti.StructuredSexual(recall_prior=True)
        networks.extend([prior_partners, sexual])
    
    # Create interventions
    interventions = []
    if params['include_testing']:
        testing = sti.HIVTest(test_prob_data=0.2, start=params['start'])
        interventions.append(testing)
    
    # Create and run simulation
    sim = sti.Sim(
        start=params['start'],
        dur=params['dur'],
        n_agents=params['n_agents'],
        diseases=hiv,
        networks=networks,
        demographics=[pregnancy, death],
        interventions=interventions
    )
    
    sim.run(verbose=1/12)
    return extract_simulation_results(sim, params)
```

#### **STIsim Components Used**

**Disease Models:**
- **`sti.HIV`**: Core HIV transmission model with configurable parameters
- **Natural History**: Acute, latent, and late-stage infection phases
- **Transmission Rates**: Male-to-female, mother-to-child, male-to-male

**Network Models:**
- **`sti.StructuredSexual`**: Structured sexual network with partner recall
- **`sti.AgeMatchedMSM`**: Age-matched men who have sex with men network
- **`sti.PriorPartners`**: Partner history tracking for recall functionality

**Interventions:**
- **`sti.HIVTest`**: HIV testing with configurable coverage
- **`sti.ART`**: Antiretroviral therapy with time-varying coverage
- **`sti.VMMC`**: Voluntary medical male circumcision
- **`sti.Prep`**: Pre-exposure prophylaxis

**Simulation Framework:**
- **`sti.Sim`**: Main simulation orchestrator
- **Demographics**: Pregnancy and death processes
- **Time Management**: Configurable start year and duration

## Data Flow Architecture

### 1. User Input Flow
```
User Interface (Sliders) → R Server → Parameter Validation → Python Function
```

### 2. Simulation Execution Flow
```
Python Function → STIsim Components → Simulation Run → Results Extraction
```

### 3. Results Processing Flow
```
Python Results → R Processing → Plotly Visualizations → User Interface
```

## File Structure

```
shiny/
├── app.R                          # Main application entry point
├── ui.R                           # User interface definition
├── server.R                       # Server logic and reactive functions
├── README.md                      # User documentation
├── ARCHITECTURE.md                # This document
├── DEPENDENCIES.md                # Dependency documentation
├── install_dependencies.R         # R package installer
├── install_all_dependencies.R     # Comprehensive installer
└── python/
    ├── hiv_simulator.py          # Python simulation engine
    └── requirements.txt          # Python dependencies
```

## Key Integration Points

### 1. R-Python Bridge
- **Entry Point**: `app.R` line 18: `source_python("python/hiv_simulator.py")`
- **Function Call**: `server.R` line 135: `py$run_hiv_simulation(params)`
- **Data Conversion**: Automatic via reticulate

### 2. Parameter Mapping
R slider inputs are mapped to STIsim parameters:
```r
# R Input → Python Parameter
input$n_agents → params['n_agents']
input$init_prev → params['init_prev']
input$beta_m2f → params['beta_m2f']
# ... etc
```

### 3. Result Processing
Python simulation results are processed for R visualization:
```python
# Python returns dictionary
return {
    'timevec': timevec,
    'years': years,
    'hiv_prevalence': hiv_prev,
    'hiv_incidence': hiv_inc,
    'population_summary': {...},
    # ... more results
}
```

```r
# R processes results for plotting
plot_ly(x = data$years, y = data$hiv_prevalence * 100, ...)
```

## Error Handling

### 1. Python Integration Errors
- **Fallback Mode**: If Python/STIsim unavailable, uses mock data
- **Error Messages**: Clear notifications to user about simulation status
- **Graceful Degradation**: App continues to function with limited features

### 2. Parameter Validation
- **Range Checking**: Ensures parameters are within valid ranges
- **Type Validation**: Converts inputs to appropriate data types
- **Error Notifications**: User-friendly error messages

### 3. Simulation Errors
- **Exception Handling**: Catches and reports simulation failures
- **Progress Tracking**: Shows simulation progress and status
- **Recovery**: Allows retry with different parameters

## Performance Considerations

### 1. Simulation Scaling
- **Population Size**: Configurable from 100 to 50,000 agents
- **Time Duration**: Up to 50 years of simulation
- **Progress Tracking**: Real-time progress updates

### 2. Memory Management
- **Efficient Data Structures**: Uses appropriate Python data types
- **Result Caching**: Stores simulation results for quick access
- **Cleanup**: Proper memory management for large simulations

### 3. User Experience
- **Responsive UI**: Non-blocking interface during simulation
- **Progress Indicators**: Visual feedback during long simulations
- **Interactive Plots**: Real-time visualization updates

## Security and Reliability

### 1. Input Sanitization
- **Parameter Validation**: All inputs validated before use
- **Range Limits**: Prevents invalid parameter values
- **Type Safety**: Ensures correct data types

### 2. Error Recovery
- **Graceful Failures**: App continues running after errors
- **User Notifications**: Clear error messages and status updates
- **Fallback Modes**: Mock data when simulation unavailable

### 3. Resource Management
- **Process Isolation**: Python simulations run in separate processes
- **Memory Limits**: Reasonable limits on simulation size
- **Timeout Handling**: Prevents infinite loops

## Deployment Considerations

### 1. Dependencies
- **R Packages**: All listed in `install_dependencies.R`
- **Python Packages**: All listed in `python/requirements.txt`
- **STIsim Package**: Installed in editable mode

### 2. Environment Setup
- **Python Virtual Environment**: Isolated Python environment
- **R Environment**: Standard R package installation
- **Cross-Platform**: Works on Windows, macOS, and Linux

### 3. Configuration
- **Port Configuration**: Configurable server port (default: 3030)
- **Host Binding**: Configurable host address (default: 0.0.0.0)
- **Logging**: Comprehensive logging for debugging

## Future Enhancements

### 1. Additional Features
- **More Disease Models**: Support for other STIs
- **Advanced Visualizations**: 3D plots, network visualizations
- **Parameter Sensitivity**: Built-in sensitivity analysis
- **Batch Processing**: Multiple simulation runs

### 2. Performance Improvements
- **Parallel Processing**: Multi-core simulation support
- **Caching**: Result caching for repeated simulations
- **Optimization**: Faster simulation algorithms

### 3. User Experience
- **Tutorial Mode**: Guided walkthrough for new users
- **Export Features**: Save results and plots
- **Collaboration**: Share simulation configurations

## Conclusion

The STIsim HIV Shiny Web App represents a sophisticated integration of modern web technologies with advanced epidemiological modeling. The architecture provides a clean separation of concerns while maintaining seamless integration between R and Python components. The application successfully bridges the gap between complex scientific modeling and user-friendly web interfaces, making advanced HIV transmission modeling accessible to a broader audience.

The modular design allows for easy maintenance, extension, and customization, while the robust error handling and fallback mechanisms ensure reliable operation across different environments and use cases.
