# STIsim Examples

Pre-configured location-specific simulations for STIsim.

This package provides `make_sim()` factory functions to create STIsim simulations with location-specific data (demographics, disease parameters, behavioral data), making it easy to get started with realistic models.

## Quick Start

```python
# Option 1: Use a location's make_sim() directly
from stisim_examples.zimbabwe.sim import make_sim
sim = make_sim()
sim.run()
sim.plot()

# Option 2: Use hivsim_examples convenience wrapper
import hivsim_examples as hx
sim = hx.Sim(location='zimbabwe')
sim.run()
sim.plot()
```

## Available Locations

- **demo** - Generic/minimal location for testing and tutorials
- **zimbabwe** - Fully parameterized Zimbabwe HIV model

List available locations programmatically:
```python
import stisim_examples as stx
print(stx.list_locations())
```

## Usage Patterns

### Pattern 1: Location make_sim() (Recommended)

Each location has a `make_sim()` function that creates a fully configured simulation:

```python
from stisim_examples.zimbabwe.sim import make_sim

sim = make_sim(n_agents=5000, dur=40)
sim.run()
```

### Pattern 2: Via hivsim_examples

```python
import hivsim_examples as hx

sim = hx.Sim(location='zimbabwe')
sim.run()
```

### Pattern 3: Base stisim with demographics + data_path

```python
import stisim as sti

sim = sti.Sim(
    location='zimbabwe',
    diseases='hiv',
    data_path='path/to/stisim_examples/zimbabwe',
    n_agents=5000,
    dur=40,
)
sim.run()
```

### Pattern 4: Base stisim with demographics only

```python
import stisim as sti

# Load only demographics, manually configure diseases
sim = sti.Sim(demographics='zimbabwe', diseases='hiv', n_agents=5000, dur=40)
sim.run()
```

### Pattern 5: Base hivsim with demographics only

```python
import hivsim as hs

# Load demographics, use default HIV configuration
sim = hs.Sim(location='zimbabwe', n_agents=5000, dur=40)
sim.run()
```

### Pattern 6: Custom data folder

```python
import stisim as sti

# Use data from external repository or custom location
sim = sti.Sim(
    location='zimbabwe',
    datafolder='/path/to/custom/data',
    diseases='hiv',
)
sim.run()
```

## What Data is Loaded?

When you use `make_sim()` for a location, the following data is loaded:

**Demographic Data** (from `stisim/data/files/`):
- Age distributions
- Fertility rates
- Death rates
- Migration patterns (if available)

**Disease-Specific Data** (from `stisim_examples/{location}/`):
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

4. **Create configuration**: `stisim_examples/southafrica/sim.py` with a `make_sim()` function

5. **Register location**: Add to `LOCATIONS` dict in `stisim_examples/loaders.py`

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

For demographic file formats, see examples in `stisim/data/files/`.

## API Reference

### `make_sim(**kwargs)` (per location)

Factory function in each location's `sim.py` that creates a fully configured simulation.

**Parameters:**
- `**kwargs`: Parameters passed to `stisim.Sim` (e.g. `n_agents`, `dur`, `interventions`)

**Returns:** Configured `stisim.Sim` instance

### `hx.Sim(location, **kwargs)` (hivsim_examples)

Convenience wrapper that calls the location's `make_sim()`.

**Parameters:**
- `location` (str): Location name ('demo', 'zimbabwe')
- `**kwargs`: Additional parameters passed to `make_sim()`

**Returns:** Configured `stisim.Sim` instance

### `stx.list_locations()`

List available locations.

**Returns:** Dictionary mapping location names to descriptions

## Testing

Run tests:
```bash
pytest tests/test_examples.py -v
```

Test specific location:
```bash
pytest tests/test_examples.py -k zim -v
```

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
