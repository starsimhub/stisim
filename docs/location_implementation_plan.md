# STIsim/HIVsim Location-Based Data Reorganization Plan

## Context

**Why**: The stisim repository needs reorganization to make it easier to create minimal examples with location-specific data. This benefits both HIV-only users (via hivsim) and multi-disease users (via stisim).

**Problem**:
- The `location` parameter in `hivsim.Sim.__init__` raises NotImplementedError (line 23 of hivsim/sim.py)
- No clear structure for location-specific data (Zimbabwe, Zambia, Kenya examples)
- stisim.Sim doesn't support location-based initialization
- Users need location data for both HIV-only and multi-disease modeling

**Goal**: Implement shared location-based data loading that works for both usage patterns:
1. `import hivsim as hsim` → `hsim.Sim(location='zambia')` - HIV only
2. `import stisim as sti` → `sti.Sim(location='zambia', diseases=['hiv', 'syph'])` - Multi-disease

Both should load the same demographic and disease-related data for consistency.

## Design Decisions

### 1. Location Data Placement: Shared in stisim/data/locations/
- Create `stisim/data/locations/` directory for all location-specific data
- Both hivsim.Sim and sti.Sim access the same location data
- Avoids duplication and ensures consistency across usage patterns
- Follows existing stisim data organization pattern (stisim/data/ already exists)

### 2. Location Data Content: Placeholder READMEs Only
- Create empty location directories (zimbabwe/, zambia/, kenya/) with comprehensive README.md files
- Include only `generic/` location with built-in Python-based parameters for testing
- Users can copy real data from external repos when needed
- Avoids data licensing questions and keeps initial implementation lean

### 3. Full Implementation of Location Loading
- Implement location loading in stisim (either extend `stisim/data/loaders.py` or create `stisim/locations.py`)
- Add `location` parameter support to `stisim.Sim.__init__`
- Make `location` parameter in `hivsim.Sim` functional (remove NotImplementedError)
- hivsim.Sim delegates to stisim's location loading
- Enable both: `hivsim.Sim(location='generic')` and `sti.Sim(location='generic', diseases=['hiv'])`

### 4. Import Patterns: Both at Root Level
- Keep hivsim at root level: `import hivsim as hsim`
- Keep stisim at root level: `import stisim as sti`
- Both packages are siblings, hivsim is a convenience wrapper for HIV-specific work
- No need to update existing test imports

### 5. Generic Parameters as Python Module
- Create `stisim/data/locations/generic/parameters.py` with functions returning DataFrames
- More flexible than CSV for minimal defaults
- Includes inline documentation
- Provides: init_prev, condom_use, art_coverage, vmmc_coverage, demographic defaults

## Target Structure

```
stisim/
├── __init__.py
├── sim.py (MODIFY - add location parameter support)
├── locations.py (NEW - ~250 lines)
├── data/
│   ├── loaders.py (existing)
│   └── locations/ (NEW)
│       ├── __init__.py (NEW - empty)
│       ├── generic/
│       │   ├── __init__.py (NEW - empty)
│       │   └── parameters.py (NEW - ~150 lines)
│       ├── zimbabwe/
│       │   └── README.md (NEW - documentation)
│       ├── zambia/
│       │   └── README.md (NEW - documentation)
│       └── kenya/
│           └── README.md (NEW - documentation)

hivsim/ (at root level, unchanged location)
├── __init__.py (MODIFY - update exports)
├── sim.py (MODIFY - use stisim's location loading)
└── README.md (UPDATE - document location usage)
```

## Implementation Steps

### Step 1: Create Location Data Directory Structure
Create the new locations directory under stisim/data/:
```bash
mkdir -p stisim/data/locations/generic
mkdir -p stisim/data/locations/zimbabwe
mkdir -p stisim/data/locations/zambia
mkdir -p stisim/data/locations/kenya
```

### Step 2: Create Generic Location Parameters
Create `/stisim/data/locations/generic/parameters.py`:
- Function `get_init_prev_hiv()` - Returns DataFrame with HIV prevalence by risk group/sex/sw
- Function `get_init_prev_syph()` - Returns DataFrame with syphilis prevalence
- Function `get_condom_use()` - Returns DataFrame with condom use by partnership type over time
- Function `get_art_coverage()` - Returns DataFrame with ART scale-up over time
- Function `get_vmmc_coverage()` - Returns DataFrame with VMMC scale-up over time
- Function `get_demographic_params()` - Returns dict with basic demographic parameters
- All functions return minimal working defaults suitable for testing

Pattern: Follow structure discovered in external repos (hiv_zim, stisim_demo_zambia, etc.):
- init_prev columns: risk_group (0-2), sex (female/male), sw (0/1), init_prev (0-1 scale)
- condom_use columns: partnership, [years as column headers], values 0-1
- coverage files: year, coverage (0-1 scale)
- Include docstrings explaining each function

### Step 3: Create Empty __init__ Files
Create placeholder files:
- `/stisim/data/locations/__init__.py` (empty)
- `/stisim/data/locations/generic/__init__.py` (empty)

### Step 4: Create Location Loader in stisim
Create `/stisim/locations.py`:
- `LOCATIONS` dict - Registry of available locations
- `Location` class - Container for location-specific data with methods:
  - `__init__(name, data_folder=None)` - Initialize and load data
  - `_load_init_prev(disease)` - Load initial prevalence for a disease (generic: Python function, others: CSV)
  - `_load_condom_use()` - Load condom use data
  - `_load_art_coverage()` - Load ART coverage
  - `_load_vmmc_coverage()` - Load VMMC coverage
  - `_load_demographic_data()` - Load demographics (age dist, fertility, deaths)
  - `get_disease_init_prev(disease_name)` - Get initial prevalence for specific disease
  - `apply_to_sim_pars(sim_pars, disease_list)` - Merge location data into sim parameters
- `load_location(name, data_folder=None)` - Factory function to create Location
- `list_locations()` - Return available locations dict

Key behavior:
- For 'generic': call Python functions from `stisim/data/locations/generic/parameters.py`
- For other locations: look for CSV files in `stisim/data/locations/{name}/`, raise FileNotFoundError with helpful message if missing
- Support multiple diseases: HIV, syphilis, gonorrhea, chlamydia, etc.
- Handle optional data (condom_use, interventions) gracefully
- Return data in format compatible with stisim modules

### Step 5: Add Location Support to stisim.Sim
Modify `/stisim/sim.py` `__init__` method:
Add location parameter handling before other initialization:
```python
def __init__(self, ..., location=None, **kwargs):
    # Handle location-based initialization
    if location is not None:
        from .locations import load_location
        loc = load_location(location)

        # Apply location data to parameters
        # This needs to merge location data for demographics and diseases
        # Based on which diseases are specified
        kwargs = loc.apply_to_sim_pars(kwargs, diseases=kwargs.get('diseases', []))

    # ... rest of existing code unchanged ...
```

### Step 6: Update hivsim.Sim to Use stisim's Location Loading
Modify `/hivsim/sim.py` line 20-24:
```python
def __init__(self, pars=None, sim_pars=None, hiv_pars=None, location=None, **kwargs):

    # Handle location-based initialization
    if location is not None:
        # Use stisim's location loading
        import stisim as sti
        loc = sti.load_location(location)

        # Merge location data into parameters
        pars = sc.mergedicts(pars, kwargs)
        pars = loc.apply_to_sim_pars(pars, diseases=['hiv'])
        kwargs = pars
        pars = None

    # ... rest of existing code unchanged ...
```

Remove the NotImplementedError and replace with functional location loading that delegates to stisim.

### Step 7: Update stisim __init__.py
Modify `/stisim/__init__.py`:
- Add after line 12: `from .locations import *`
- This exports load_location, list_locations, Location at stisim package level

### Step 8: Update hivsim __init__.py
Modify `/hivsim/__init__.py`:
- Keep existing: `from .sim import *`
- Add: `import stisim; load_location = stisim.load_location; list_locations = stisim.list_locations`
- Or add to __all__: `'load_location', 'list_locations'`
- This makes location functions available via hivsim too

### Step 9: Create Location Documentation
Create comprehensive README.md files for each location:

**`/stisim/data/locations/zimbabwe/README.md`**:
- Document required CSV file formats (init_prev_hiv.csv, condom_use.csv, etc.)
- Show example rows for each file type
- List data sources and references
- Point to external repo: `/Users/robynstuart/gf/hiv_zim/data/`

Repeat for zambia/ and kenya/ with appropriate external repo paths.

Include example file formats and point to external data sources.

### Step 10: Add Location Loading Test
Add to `/tests/test_hiv.py`:
```python
def test_location_loading_hivsim():
    """ Test location loading functionality via hivsim """
    import hivsim

    # Test generic location
    loc = hivsim.load_location('generic')
    assert loc.get_disease_init_prev('hiv') is not None
    assert loc.condom_data is not None

    # Test creating sim with location
    sim = hivsim.Sim(location='generic', n_agents=100, dur=10)
    sim.run()
    assert 'hiv' in sim.diseases

    # Test listing locations
    locs = hivsim.list_locations()
    assert 'generic' in locs
    assert 'zimbabwe' in locs

    return sim


def test_location_loading_stisim():
    """ Test location loading functionality via stisim """
    import stisim as sti

    # Test generic location
    loc = sti.load_location('generic')
    assert loc.get_disease_init_prev('hiv') is not None

    # Test creating HIV sim with location
    sim1 = sti.Sim(location='generic', diseases=['hiv'], n_agents=100, dur=10)
    sim1.run()
    assert 'hiv' in sim1.diseases

    # Test creating multi-disease sim with location
    sim2 = sti.Sim(location='generic', diseases=['hiv', 'syphilis'], n_agents=100, dur=10)
    sim2.run()
    assert 'hiv' in sim2.diseases
    assert 'syphilis' in sim2.diseases

    return sim1, sim2
```

### Step 11: Update Documentation
**`/hivsim/README.md`**:
- Add "Quick Start" section with import examples
- Add "Location Data" section explaining available locations and usage
- Add code examples for location-based sims
- Explain relationship to stisim

**`/README.md`**:
- Add brief section on location-based simulations
- Show examples for both hivsim and stisim usage patterns
- Reference location data documentation in `stisim/data/locations/`

## Critical Files

### New Files (7 files)
1. `/stisim/locations.py` - Core location loading logic (~250 lines)
2. `/stisim/data/locations/generic/parameters.py` - Built-in test parameters (~150 lines)
3. `/stisim/data/locations/__init__.py` - Empty
4. `/stisim/data/locations/generic/__init__.py` - Empty
5. `/stisim/data/locations/zimbabwe/README.md` - Documentation (~70 lines)
6. `/stisim/data/locations/zambia/README.md` - Documentation (~70 lines)
7. `/stisim/data/locations/kenya/README.md` - Documentation (~70 lines)

### Modified Files (6 files)
1. `/stisim/__init__.py` - Add `from .locations import *`
2. `/stisim/sim.py` - Add location parameter support
3. `/hivsim/__init__.py` - Add location function exports/aliases
4. `/hivsim/sim.py` - Implement location parameter (replace NotImplementedError, delegate to stisim)
5. `/hivsim/README.md` - Update with location usage examples
6. `/tests/test_hiv.py` - Add location loading tests

### No Deleted Files
- hivsim/ stays at root level
- No file moves required

## Verification Steps

### 1. Test Import Syntax for Both Patterns
```python
# Pattern 1: HIV-only via hivsim
import hivsim as hsim
sim1 = hsim.Sim(location='generic', n_agents=100, dur=5)
sim1.run()

# Pattern 2: Multi-disease via stisim
import stisim as sti
sim2 = sti.Sim(location='generic', diseases=['hiv'], n_agents=100, dur=5)
sim2.run()

# Pattern 3: Multi-disease HIV+Syphilis
sim3 = sti.Sim(location='generic', diseases=['hiv', 'syphilis'], n_agents=100, dur=5)
sim3.run()
```

### 2. Run Existing Tests
```bash
cd tests
pytest test_hiv.py -v
pytest devtests/devtest_hiv.py -v
```
All existing tests should pass unchanged (no import updates needed).

### 3. Test Location Loading via hivsim
```python
import hivsim

# List locations
print(hivsim.list_locations())

# Load generic location
loc = hivsim.load_location('generic')
print("HIV init prev:", loc.get_disease_init_prev('hiv').head())
print("Condom use:", loc.condom_data.head())

# Create sim with location
sim = hivsim.Sim(location='generic', n_agents=100, dur=10)
sim.run()
print("HIV prevalence:", sim.results.hiv.prevalence[-1])
```

### 4. Test Location Loading via stisim
```python
import stisim as sti

# Load location
loc = sti.load_location('generic')
print("HIV init prev:", loc.get_disease_init_prev('hiv').head())

# HIV-only sim
sim1 = sti.Sim(location='generic', diseases=['hiv'], n_agents=100, dur=10)
sim1.run()

# Multi-disease sim
sim2 = sti.Sim(location='generic', diseases=['hiv', 'syphilis'], n_agents=100, dur=10)
sim2.run()
print("HIV prevalence:", sim2.results.hiv.prevalence[-1])
print("Syphilis prevalence:", sim2.results.syphilis.prevalence[-1])
```

### 5. Test Without Location (Should Still Work)
```python
import hivsim
import stisim as sti

# These should work without location parameter (existing functionality)
sim1 = hivsim.Sim(n_agents=100, dur=5)
sim1.run()

sim2 = sti.Sim(diseases=['hiv'], n_agents=100, dur=5)
sim2.run()
```

### 6. Run Full Test Suite
```bash
cd tests
pytest test_*.py -n auto
```
Ensure no regressions in other parts of stisim.

### 7. Verify Package Installation
```bash
pip install -e .
python -c "import hivsim; print(hivsim.__file__)"
python -c "import hivsim; print(hivsim.list_locations())"
python -c "import stisim; print(stisim.list_locations())"
```

### 8. Check Documentation
- Read through updated READMEs
- Verify code examples work as written
- Check that location READMEs are clear and complete

## Success Criteria

1. ✅ All existing tests pass with NO changes required
2. ✅ Both import patterns work: `import hivsim as hsim` and `import stisim as sti`
3. ✅ Location loading works for both: `hsim.Sim(location='generic')` and `sti.Sim(location='generic', diseases=['hiv'])`
4. ✅ Multi-disease location loading works: `sti.Sim(location='generic', diseases=['hiv', 'syphilis'])`
5. ✅ Sims work without location parameter (existing functionality preserved)
6. ✅ New tests `test_location_loading_hivsim()` and `test_location_loading_stisim()` pass
7. ✅ Documentation is complete and examples work
8. ✅ No regression in other stisim functionality
9. ✅ Package installs cleanly with `pip install -e .`

## Implementation Notes

### Data Format Reference
Based on exploration of 6 external repos, location data follows these patterns:

**init_prev_hiv.csv**:
```csv
risk_group,sex,sw,init_prev
0,female,0,0.001
1,female,0,0.01
2,female,0,0.05
0,female,1,0.15
...
```

**condom_use.csv**:
```csv
partnership,1970,1980,1990,2000,2010,2020
"(0,0)",0,0,0,0.05,0.1,0.15
"(1,1)",0,0,0.1,0.5,0.7,0.8
"(fsw,client)",0,0.1,0.5,0.7,0.85,0.9
...
```

**n_art.csv** / **n_vmmc.csv**:
```csv
year,coverage
2000,0
2005,0.1
2010,0.4
2020,0.8
```

### Key Dependencies
- Pattern follows existing `stisim/data/loaders.py`
- HIV parameters defined in `stisim/diseases/hiv.py` (HIVPars class)
- Network parameters for StructuredSexual in `stisim/networks.py`
- Location data discovered in external repos: hiv_kenya, hiv_zambia, hiv_zim, stisim_vddx_zim, stisim_demo_zambia, syph_dx_zim

### Risk Mitigation
- **Risk**: Location loading errors with missing data
  - **Mitigation**: Clear error messages pointing to README.md for file format specs
- **Risk**: Confusion about shared location data between hivsim and stisim
  - **Mitigation**: Comprehensive documentation explaining that location data is stored in stisim/data/locations/ and shared between both
- **Risk**: Conflicts between location parameters and manual disease/demographic parameters
  - **Mitigation**: Location data should be merged carefully, allowing user overrides
- **Risk**: Existing tests may break if they don't handle new location parameter
  - **Mitigation**: Location parameter is optional, defaults to None, so existing code continues to work

## Future Enhancements

1. Populate zimbabwe/, zambia/, kenya/ with real data from external repos
2. Add data validation utilities to check CSV format correctness
3. Integrate with data downloaders for automatic data fetching
4. Add more locations (South Africa, Uganda, etc.)
5. Create calibration workflow for new locations
6. Add visualization tools for location data exploration
