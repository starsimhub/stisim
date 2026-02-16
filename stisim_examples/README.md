# STIsim Examples

Pre-configured location-specific simulations for STIsim.

This package provides factory functions to create STIsim simulations with location-specific data (demographics, disease parameters, behavioral data), making it easy to get started with realistic models.

## Quick Start

```python
import stisim_examples as stx

# Zimbabwe HIV model
sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
sim.run()
sim.plot()

# Or use the HIV-specific convenience function
sim = stx.HIVSim(location='zimbabwe')
sim.run()
```

## Available Locations

- **demo** - Generic/minimal location for testing and tutorials
- **zimbabwe** - Fully parameterized Zimbabwe HIV model
- **kenya** - Placeholder (data not yet populated)

List available locations programmatically:
```python
print(stx.list_locations())
```

## Usage Patterns

### Pattern 1: Full Pre-configured Examples (Recommended)

```python
import stisim_examples as stx

# Single disease
sim = stx.Sim(demographics='zimbabwe', diseases='hiv', n_agents=5000, dur=40)
sim.run()

# Multiple diseases (when data is available)
sim = stx.Sim(demographics='zimbabwe', diseases=['hiv', 'syphilis'])
sim.run()
```

### Pattern 2: HIV-Specific Convenience

```python
import stisim_examples as stx

# Cleaner syntax for HIV-only modeling
sim = stx.HIVSim(location='zimbabwe', n_agents=5000, dur=40)
sim.run()
```

### Pattern 3: Via hivsim_examples

```python
import hivsim_examples as hx

# Same as stx.HIVSim
sim = hx.Sim(location='zimbabwe')
sim.run()
```

### Pattern 4: Base stisim with Demographics Only

```python
import stisim as sti

# Load only demographics, manually configure diseases
sim = sti.Sim(demographics='zimbabwe', diseases='hiv', n_agents=5000, dur=40)
sim.run()
```

### Pattern 5: Base hivsim with Demographics Only

```python
import hivsim as hs

# Load demographics, use default HIV configuration
sim = hs.Sim(location='zimbabwe', n_agents=5000, dur=40)
sim.run()
```

### Pattern 6: Custom Data Folder

```python
import stisim as sti

# Use data from external repository or custom location
sim = sti.Sim(
    demographics='zimbabwe',
    datafolder='/Users/username/custom/data',
    diseases='hiv'
)
sim.run()
```

## What Data is Loaded?

When you use `stx.Sim(demographics='zimbabwe', diseases='hiv')`, the following data is loaded:

**Demographic Data** (from `stisim/data/files/`):
- Age distributions
- Fertility rates
- Death rates
- Migration patterns (if available)

**Disease-Specific Data** (from `stisim_examples/zimbabwe/`):
- Initial HIV prevalence by risk group/sex
- Condom use by partnership type over time
- ART coverage scale-up
- VMMC coverage scale-up

Disease-specific data is fully integrated into the simulation:
- **Initial prevalence** data is passed to `sti.HIV(init_prev_data=...)` for risk-group/sex/SW-stratified initialization
- **Condom use** data is passed to `sti.StructuredSexual(condom_data=...)` for time-varying condom use by partnership type
- **ART coverage** data is passed to `sti.ART(coverage_data=...)` for data-driven ART scale-up
- **VMMC coverage** data is passed to `sti.VMMC(coverage_data=...)` for data-driven VMMC scale-up
- **HIV testing** interventions are created with default scale-up curves for FSW, general population, and low-CD4 testing

## Adding New Locations

To add a new location (e.g., "southafrica"):

1. **Create directory**: `stisim_examples/southafrica/`

2. **Add STI data files**:
   - `init_prev_hiv.csv` - HIV prevalence by risk group/sex/sw
   - `condom_use.csv` - Condom use by partnership and year
   - `art_coverage.csv` - ART coverage over time
   - `vmmc_coverage.csv` - VMMC coverage over time

3. **Add demographic data** to `stisim/data/files/`:
   - `southafrica_asfr.csv` - Age-specific fertility rates
   - `southafrica_deaths.csv` - Death rates
   - `southafrica_age_YYYY.csv` - Age distribution

4. **Create configuration**: `stisim_examples/southafrica/sim.py`

5. **Register location**: Add to `LOCATIONS` dict in `stisim_examples/loaders.py`

6. **Document**: Create `stisim_examples/southafrica/README.md`

See `zimbabwe/` for a complete example, and `demo/` for minimal file formats.

## File Format Specifications

### init_prev_hiv.csv
```csv
risk_group,sex,sw,init_prev
0,female,0,0.001
1,female,0,0.01
2,female,0,0.05
0,female,1,0.15
...
```

### condom_use.csv
```csv
partnership,1980,1990,2000,2010,2020
"(0,0)",0.00,0.00,0.05,0.10,0.15
"(1,1)",0.00,0.10,0.30,0.50,0.60
"(fsw,client)",0.05,0.30,0.60,0.80,0.90
...
```

### art_coverage.csv / vmmc_coverage.csv
```csv
year,coverage
2000,0.00
2005,0.10
2010,0.40
...
```

For demographic file formats, see examples in `stisim/data/files/` or `tests/test_data/`.

## API Reference

### `stx.Sim(demographics, diseases, **kwargs)`

Factory function to create pre-configured simulations.

**Parameters:**
- `demographics` (str): Location name ('demo', 'zimbabwe', 'kenya')
- `diseases` (str/list): Disease(s) to include ('hiv', ['hiv', 'syphilis'])
- `**kwargs`: Additional parameters passed to `stisim.Sim`

**Returns:** Configured `stisim.Sim` instance

### `stx.HIVSim(location, **kwargs)`

Convenience function for HIV-only simulations.

**Parameters:**
- `location` (str): Location name
- `**kwargs`: Additional parameters

**Returns:** Configured `stisim.Sim` instance with HIV

### `stx.list_locations()`

List available locations.

**Returns:** Dictionary mapping location names to descriptions

## Testing

Run tests:
```bash
cd tests
pytest test_stisim_examples.py -v
```

Test specific location:
```bash
pytest test_stisim_examples.py -k zimbabwe -v
```

## Contributing

Contributions of new location data are welcome! Please ensure:
- Data files follow format specifications
- Configuration helpers are provided in `{location}/sim.py`
- Comprehensive README documents data sources
- Tests are added to verify functionality

## Related Packages

- **stisim** - Core STI modeling framework
- **hivsim** - HIV-specific modeling convenience package
- **starsim** - Underlying agent-based modeling framework

## License

Same as STIsim (MIT License)

## Support

For questions or issues:
- GitHub: https://github.com/starsimhub/stisim/issues
- Documentation: https://docs.idmod.org/projects/stisim
