# STIsim Examples Implementation Plan

## Overview

Create a `stisim_examples` package alongside `stisim` that provides location-specific, pre-configured simulations. Uses Starsim's `register_modules()` pattern to enable string-based initialization.

## Design

### Key Principles
1. **Factory function pattern**: `stx.Sim()` is a factory that returns configured `stisim.Sim` instances
2. **Location-based data**: Each location folder contains CSV data files and helper functions
3. **Shared demographics**: Demographic data lives in `stisim/data/`, disease/behavioral data in `stisim_examples/{location}/`
4. **Dual usage**: Full examples via `stx.Sim()` OR base stisim with just demographics
5. **Convenience wrappers**: `stx.HIVSim()` and `hivsim_examples` for HIV-only workflows

## Target Structure

```
stisim/
├── data/
│   ├── demographics/       # Shared demographic data
│   └── ...

stisim_examples/            # NEW - Top-level package
├── __init__.py             # NEW - Factory functions (Sim, HIVSim)
├── loaders.py              # NEW - Data loading utilities
├── zimbabwe/
│   ├── __init__.py         # NEW
│   ├── sim.py              # NEW - Helper functions
│   ├── init_prev_hiv.csv   # NEW
│   ├── init_prev_syph.csv  # NEW
│   ├── condom_use.csv      # NEW
│   ├── art_coverage.csv    # NEW
│   └── README.md           # NEW - Documentation
├── kenya/
│   ├── __init__.py
│   ├── sim.py
│   └── README.md           # Data files TBD
└── demo/
    ├── __init__.py
    ├── sim.py
    ├── *.csv
    └── README.md

hivsim_examples/            # NEW - Thin wrapper package
├── __init__.py             # NEW - Aliases stx.HIVSim

tests/
├── test_hiv.py             # MODIFY - Add stisim_examples tests
└── test_stisim_examples.py # NEW - Comprehensive examples tests
```

## Usage Patterns

### Pattern 1: Full Pre-configured Examples
```python
import stisim_examples as stx
import starsim as ss

# Register so we can use string references (optional but useful)
ss.register_modules(stx)

# Zimbabwe HIV model
sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
sim.run()

# Zimbabwe HIV-Syphilis model
sim = stx.Sim(demographics='zimbabwe', diseases=['hiv', 'syphilis'])
sim.run()

# Kenya HIV model
sim = stx.Sim(demographics='kenya', diseases='hiv')
sim.run()
```

### Pattern 2: HIV-specific convenience
```python
import stisim_examples as stx

# Cleaner syntax for HIV-only
sim = stx.HIVSim(location='zimbabwe')
sim.run()
```

### Pattern 3: hivsim_examples alias
```python
import hivsim_examples as hx

# Same as stx.HIVSim
sim = hx.Sim(location='zimbabwe')
sim.run()
```

### Pattern 4: Base stisim with demographics from custom path
```python
import stisim as sti

# Load demographics from custom data folder
sim = sti.Sim(
    demographics='zimbabwe',
    data_folder='/Users/robynstuart/gf/hiv_zim/data',
    diseases='hiv'
)
sim.run()

# Or use default stisim/data/ location
sim = sti.Sim(demographics='zimbabwe', diseases='hiv')
sim.run()
```

### Pattern 5: Base hivsim with demographics from custom path
```python
import hivsim as hs

# Load demographics from custom data folder
sim = hs.Sim(
    demographics='zimbabwe',
    data_folder='/Users/robynstuart/gf/hiv_zim/data'
)
sim.run()

# Or use default stisim/data/ location
sim = hs.Sim(demographics='zimbabwe')
sim.run()
```

## Implementation Steps

### Step 1: Create stisim_examples directory structure

```bash
mkdir -p stisim_examples/zimbabwe
mkdir -p stisim_examples/kenya
mkdir -p stisim_examples/demo
```

### Step 2: Create stisim_examples/loaders.py

Core data loading utilities (~150 lines):

```python
import pandas as pd
import sciris as sc
from pathlib import Path

__all__ = ['load_location_data', 'LOCATIONS']

# Registry of available locations
LOCATIONS = {
    'zimbabwe': 'Zimbabwe',
    'kenya': 'Kenya',
    'demo': 'Demo/Generic location',
}

def load_location_data(location, diseases=None):
    """
    Load location-specific data for specified diseases.

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        diseases (str/list): Disease(s) to load data for

    Returns:
        dict: Nested dict with disease parameters and behavioral data
    """
    # Normalize diseases to list
    if isinstance(diseases, str):
        diseases = [diseases]

    # Get location path
    loc_path = Path(__file__).parent / location
    if not loc_path.exists():
        raise ValueError(f"Location '{location}' not found")

    # Load data
    data = sc.objdict()
    data.location = location
    data.diseases = {}

    # Load disease-specific data
    for disease in diseases:
        disease_data = sc.objdict()

        # Load initial prevalence
        init_prev_file = loc_path / f'init_prev_{disease}.csv'
        if init_prev_file.exists():
            disease_data.init_prev = pd.read_csv(init_prev_file)

        # Add disease-specific loaders here
        # e.g., ART coverage for HIV, treatment for syphilis, etc.

        data.diseases[disease] = disease_data

    # Load behavioral data (shared across diseases)
    condom_file = loc_path / 'condom_use.csv'
    if condom_file.exists():
        data.condom_use = pd.read_csv(condom_file)

    # Load intervention data
    art_file = loc_path / 'art_coverage.csv'
    if art_file.exists():
        data.art_coverage = pd.read_csv(art_file)

    return data
```

### Step 3: Create stisim_examples/__init__.py

Factory functions for Sim creation (~100 lines):

```python
import stisim as sti
from .loaders import load_location_data, LOCATIONS

__all__ = ['Sim', 'HIVSim', 'list_locations']

def list_locations():
    """List available locations."""
    return LOCATIONS

def Sim(demographics=None, diseases=None, **kwargs):
    """
    Factory function to create pre-configured stisim simulations.

    Args:
        demographics (str): Location name for demographics (zimbabwe, kenya, demo)
        diseases (str/list): Disease(s) to include
        **kwargs: Additional parameters passed to stisim.Sim

    Returns:
        stisim.Sim: Configured simulation instance

    Examples:
        >>> sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
        >>> sim = stx.Sim(demographics='kenya', diseases=['hiv', 'syphilis'])
    """
    if demographics is None:
        raise ValueError("demographics parameter is required")

    if diseases is None:
        raise ValueError("diseases parameter is required")

    # Load location-specific data
    loc_data = load_location_data(demographics, diseases)

    # Import location-specific helper functions
    # e.g., from .zimbabwe.sim import configure_zimbabwe_hiv
    # and use them to configure parameters

    # Merge location data into kwargs
    # This is where the magic happens - converting CSV data into
    # disease parameters, network parameters, etc.

    # For now, placeholder - detailed implementation in Step 6
    pars = kwargs.copy()
    pars['demographics'] = demographics
    pars['diseases'] = diseases

    # Create and return the sim
    sim = sti.Sim(**pars)

    return sim

def HIVSim(location=None, **kwargs):
    """
    Convenience function for HIV-only simulations.

    Args:
        location (str): Location name (zimbabwe, kenya, demo)
        **kwargs: Additional parameters passed to Sim()

    Returns:
        stisim.Sim: Configured HIV simulation

    Example:
        >>> sim = stx.HIVSim(location='zimbabwe')
    """
    if location is None:
        raise ValueError("location parameter is required")

    return Sim(demographics=location, diseases='hiv', **kwargs)
```

### Step 4: Create location helper modules

**stisim_examples/zimbabwe/sim.py** (~100 lines):

```python
"""
Zimbabwe-specific simulation configuration helpers.
"""
import pandas as pd
import sciris as sc
from pathlib import Path

def load_zimbabwe_data(diseases):
    """Load all Zimbabwe data files."""
    from ..loaders import load_location_data
    return load_location_data('zimbabwe', diseases)

def configure_hiv_pars(data):
    """
    Configure HIV parameters from Zimbabwe data.

    Args:
        data: Data dict from load_location_data()

    Returns:
        dict: HIV parameters for stisim
    """
    hiv_pars = {}

    # Convert init_prev CSV to format needed by stisim.HIV
    if 'init_prev' in data.diseases.hiv:
        init_prev_df = data.diseases.hiv.init_prev
        # Transform DataFrame into parameter structure
        # This depends on stisim.HIV's expected format
        hiv_pars['init_prev'] = init_prev_df

    return hiv_pars

def configure_network_pars(data):
    """Configure network parameters from behavioral data."""
    net_pars = {}

    if 'condom_use' in data:
        # Transform condom use data
        net_pars['condom_use'] = data.condom_use

    return net_pars

# Similar functions for syphilis, other diseases
```

Repeat for `kenya/sim.py` and `demo/sim.py`.

### Step 5: Create hivsim_examples package

**hivsim_examples/__init__.py** (~10 lines):

```python
"""
HIV-specific examples package.
Thin wrapper around stisim_examples.HIVSim.
"""
import stisim_examples as stx

# Alias HIVSim as Sim for cleaner syntax
Sim = stx.HIVSim
list_locations = stx.list_locations

__all__ = ['Sim', 'list_locations']
```

### Step 6: Add demographics parameter to stisim.Sim

**Modify stisim/sim.py** (~30 lines of changes):

```python
def __init__(self, ..., demographics=None, data_folder=None, **kwargs):
    """
    Initialize simulation.

    Args:
        demographics (str): Optional location name to load demographic data
        data_folder (str/Path): Optional custom path to data folder.
                                If None, uses stisim/data/
        ...
    """
    # Handle demographics parameter
    if demographics is not None:
        # Load demographic data from specified or default location
        demo_data = self._load_demographics(demographics, data_folder)
        kwargs = sc.mergedicts(demo_data, kwargs)

    # ... rest of existing __init__ code ...

def _load_demographics(self, location, data_folder=None):
    """
    Load demographic data for a location.

    Args:
        location (str): Location name (e.g., 'zimbabwe')
        data_folder (str/Path): Custom data folder, or None for default stisim/data/

    Returns:
        dict: Demographic parameters
    """
    from pathlib import Path
    from .data import loaders  # Or wherever demographic loaders live

    # Use custom path or default to stisim/data/
    if data_folder is None:
        data_folder = Path(__file__).parent / 'data'
    else:
        data_folder = Path(data_folder)

    return loaders.load_demographics(location, data_folder)
```

### Step 7: Add demographics parameter to hivsim.Sim

**Modify hivsim/sim.py** (~20 lines):

```python
def __init__(self, pars=None, sim_pars=None, hiv_pars=None, demographics=None, data_folder=None, **kwargs):
    """
    Args:
        demographics (str): Optional location name to load demographic data
        data_folder (str/Path): Optional custom path to data folder
        ...
    """
    # Pass demographics and data_folder through to stisim.Sim
    if demographics is not None:
        kwargs['demographics'] = demographics
    if data_folder is not None:
        kwargs['data_folder'] = data_folder

    # ... rest of existing __init__ code ...
```

### Step 8: Populate Zimbabwe data files

Copy data from external repos:
- From `/Users/robynstuart/gf/hiv_zim/data/`
- Create `init_prev_hiv.csv`, `condom_use.csv`, `art_coverage.csv`, etc.
- Include comprehensive README.md with data sources

**stisim_examples/zimbabwe/README.md**:

```markdown
# Zimbabwe HIV Model Data

## Data Sources

- **HIV initial prevalence**: [Source]
- **Condom use**: [Source]
- **ART coverage**: [Source]

## File Formats

### init_prev_hiv.csv
Columns: risk_group, sex, sw, init_prev
...

### condom_use.csv
Columns: partnership, [years]
...
```

### Step 9: Create demo location

**stisim_examples/demo/sim.py**: Simple generic parameters (~50 lines)

Create minimal working defaults suitable for testing/tutorials.

### Step 10: Add comprehensive tests

**tests/test_stisim_examples.py** (~200 lines):

```python
import stisim_examples as stx
import hivsim_examples as hx
import starsim as ss

def test_stx_list_locations():
    """Test listing available locations."""
    locs = stx.list_locations()
    assert 'zimbabwe' in locs
    assert 'kenya' in locs
    assert 'demo' in locs

def test_stx_sim_zimbabwe_hiv():
    """Test Zimbabwe HIV simulation."""
    sim = stx.Sim(demographics='zimbabwe', diseases='hiv', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    assert sim.results.hiv.prevalence[-1] > 0

def test_stx_sim_zimbabwe_multi():
    """Test Zimbabwe multi-disease simulation."""
    sim = stx.Sim(demographics='zimbabwe', diseases=['hiv', 'syphilis'], n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases
    assert 'syphilis' in sim.diseases

def test_stx_hivsim():
    """Test HIVSim convenience function."""
    sim = stx.HIVSim(location='zimbabwe', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases

def test_hx_sim():
    """Test hivsim_examples wrapper."""
    sim = hx.Sim(location='zimbabwe', n_agents=100, dur=5)
    sim.run()
    assert 'hiv' in sim.diseases

def test_register_modules():
    """Test register_modules pattern."""
    ss.register_modules(stx)
    # After registration, string-based init should work
    # (if we implement it that way)

def test_base_stisim_demographics():
    """Test base stisim with demographics only."""
    import stisim as sti
    sim = sti.Sim(demographics='zimbabwe', diseases='hiv', n_agents=100, dur=5)
    sim.run()

def test_base_hivsim_demographics():
    """Test base hivsim with demographics only."""
    import hivsim as hs
    sim = hs.Sim(demographics='zimbabwe', n_agents=100, dur=5)
    sim.run()

def test_custom_data_folder():
    """Test custom data folder path."""
    import stisim as sti
    # Assuming data exists at this path
    sim = sti.Sim(
        demographics='zimbabwe',
        data_folder='/Users/robynstuart/gf/hiv_zim/data',
        diseases='hiv',
        n_agents=100,
        dur=5
    )
    sim.run()

def test_hivsim_custom_data_folder():
    """Test custom data folder with hivsim."""
    import hivsim as hs
    sim = hs.Sim(
        demographics='zimbabwe',
        data_folder='/Users/robynstuart/gf/hiv_zim/data',
        n_agents=100,
        dur=5
    )
    sim.run()
```

**Modify tests/test_hiv.py** - Add subset of above tests to existing test suite.

### Step 11: Update documentation

**README.md** - Add examples section:

```markdown
## Location-Specific Examples

STIsim includes pre-configured location-specific models via `stisim_examples`:

```python
import stisim_examples as stx

# Zimbabwe HIV model with location-specific data
sim = stx.Sim(demographics='zimbabwe', diseases='hiv')
sim.run()

# Or use the HIV-specific convenience function
sim = stx.HIVSim(location='zimbabwe')
```

See `stisim_examples/` for available locations and data.
```

**stisim_examples/README.md** - Create comprehensive guide (~100 lines).

## Implementation Order

1. ✅ **Step 1**: Create directory structure
2. ✅ **Step 2**: Create loaders.py with data loading utilities
3. ✅ **Step 3**: Create stisim_examples/__init__.py with factory functions
4. ✅ **Step 9**: Create demo location (simple test data)
5. ✅ **Step 4**: Create demo/sim.py helper functions
6. ✅ **Step 6**: Add demographics to stisim.Sim
7. ✅ **Step 7**: Add demographics to hivsim.Sim
8. ✅ **Step 10**: Add tests for demo location
9. ✅ **Step 5**: Create hivsim_examples package
10. ✅ **Step 10**: Add tests for hivsim_examples
11. ✅ **Step 8**: Populate Zimbabwe data
12. ✅ **Step 4**: Create zimbabwe/sim.py helpers
13. ✅ **Step 10**: Add tests for Zimbabwe
14. ✅ **Step 8**: Populate Kenya data (or add placeholder)
15. ✅ **Step 11**: Update documentation

## Success Criteria

1. ✅ `import stisim_examples as stx` works
2. ✅ `stx.Sim(demographics='zimbabwe', diseases='hiv')` creates configured sim
3. ✅ `stx.HIVSim(location='zimbabwe')` works as shorthand
4. ✅ `import hivsim_examples as hx` and `hx.Sim(location='zimbabwe')` works
5. ✅ `sti.Sim(demographics='zimbabwe', diseases='hiv')` loads demographics from stisim/data/
6. ✅ `sti.Sim(demographics='zimbabwe', data_folder='/custom/path', diseases='hiv')` loads from custom path
7. ✅ `hs.Sim(demographics='zimbabwe')` loads demographics from stisim/data/
8. ✅ `hs.Sim(demographics='zimbabwe', data_folder='/custom/path')` loads from custom path
9. ✅ All tests pass including new stisim_examples tests
8. ✅ No regression in existing tests
9. ✅ Documentation is complete with working examples
10. ✅ Zimbabwe data is properly loaded and integrated

## Key Design Decisions

### Why factory functions instead of classes?
- Simpler implementation
- Cleaner syntax: `stx.Sim()` vs `stx.Sim()`
- Returns standard stisim.Sim instances (no subclassing needed)
- More flexible for future enhancements

### Why separate stisim_examples package?
- Clear separation between core library and examples
- Examples can have different dependencies/data
- Easier to maintain and extend
- Users can choose: minimal install (stisim only) or with examples

### Why demographics in stisim/data/ but disease data in stisim_examples/?
- Demographics are fundamental and reusable
- Disease/behavioral data is more example-specific
- Allows base stisim to support demographics without full examples
- Balances code reuse with separation of concerns

### Why hivsim_examples as alias?
- Consistency with existing hivsim package
- Clearer intent for HIV-only users
- Minimal code duplication (just imports)
- Can evolve independently if needed

## Future Enhancements

1. Add more locations (South Africa, Uganda, etc.)
2. Data validation utilities
3. Calibration workflows for new locations
4. Interactive data exploration tools
5. String-based initialization via register_modules pattern
6. Automatic data downloading/updating from external sources
7. Support custom data paths in stisim_examples: `stx.Sim(demographics='zimbabwe', data_folder='/custom/path', diseases='hiv')`
