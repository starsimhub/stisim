# HIVsim Example Projects Plan

## Context

**Why**: Make it easy to create and share location-specific HIV models for testing, demonstration, and as starting points for real projects. Dev priority is HIVsim right now.

**Two problems** (we're solving #1 now, keeping #2 in mind):

1. **Example projects** (this plan): Ship a couple of pre-configured HIV examples in stisim so users can do `hs.demo('zimbabwe')` and get a working, data-driven HIV sim. These also serve as templates for users building their own country models.

2. **Automated data loading** (future): Eventually support `hs.Sim(location='kenya')` backed by an external `stisim_data` repo with auto-extracted demographic and disease data per country. Not in scope now, but the example project structure should be forward-compatible with it.

## What Already Exists (v1.4.8)

These pieces are already in place:

- **`sti.Sim`** accepts `location` and `data_path` parameters. When `data_path` is provided, it uses `DataLoader` to create disease/network/intervention modules from CSVs.
- **`DataLoader`** (in `stisim/data/loaders.py`) loads `init_prev_*.csv`, `condom_use.csv`, `art_coverage.csv`, `vmmc_coverage.csv`, and `{location}_hiv_data.csv` from a folder, then creates configured module instances.
- **`hivsim.Sim`** passes `location` through to `sti.Sim` as `demographics`.
- **Demographic loading** in `sti.Sim.process_demographics()` loads age distributions, fertility, mortality, and migration from CSV files when `demographics` is a string.
- **Plotting** via starsim's auto_plot system: `sim.plot('hiv')` produces a curated panel with data overlay when `sim.data` is set.
- **External example repos** (hiv_kenya, hiv_zambia, hiv_zim, stisim_demo_zambia, stisim_vddx_zim) follow a consistent pattern: `data/` folder with CSVs + `make_sim()` function.

## Design

### User-facing API

```python
import hivsim as hs

# Quick demo — creates, runs, and plots a pre-configured HIV sim
hs.demo()              # Simple default HIV sim (no data files needed)
hs.demo('zimbabwe')    # Zimbabwe HIV model with real data + data overlay

# Get the sim without running it
sim = hs.demo('zimbabwe', run=False)
sim.run()
sim.plot()             # Auto-plots HIV results with UNAIDS data overlay

# List what's available
hs.demo(list=True)     # Prints available examples

# For customization, import the example's make_sim directly
from hivsim_examples.zimbabwe.sim import make_sim
sim = make_sim(n_agents=5000, dur=40)
sim.run()
```

This mirrors `ss.demo()` in starsim, but with an optional example name argument. The `demo()` function lives in hivsim (not stisim) because the disease is always HIV — stisim is disease-agnostic so `sti.demo()` would be ambiguous.

### hivsim_examples package structure

```
hivsim_examples/
├── __init__.py              # EXAMPLES registry, get_example()
├── simple/
│   ├── __init__.py
│   └── sim.py               # make_sim(**kwargs) → hs.Sim (no CSV files)
└── zimbabwe/
    ├── __init__.py
    ├── sim.py               # make_sim(**kwargs) → sti.Sim with data_path
    ├── init_prev_hiv.csv
    ├── condom_use.csv
    ├── art_coverage.csv
    ├── vmmc_coverage.csv
    └── zimbabwe_hiv_data.csv
```

### The two examples

#### `simple` — Minimal, no data files

The simple example shows the most basic way to set up an HIV sim. No CSV data files, no `DataLoader` — just `hs.Sim()` with hivsim defaults:

```python
import hivsim as hs

def make_sim(**kwargs):
    sim = hs.Sim(
        n_agents=2000,
        dur=20,
        **kwargs,
    )
    return sim
```

`hivsim.Sim` provides sensible defaults: HIV disease, StructuredSexual network, basic testing/ART/VMMC/PrEP, Pregnancy + Deaths demographics. Users can see exactly what a minimal working sim looks like.

#### `zimbabwe` — Full data-driven country model

The Zimbabwe example shows the full pattern for a real country model. Uses `data_path` + `DataLoader` to load CSVs, has calibrated parameters, custom testing interventions, and UNAIDS data for plot overlay:

```python
sim_pars = dict(start=1990, stop=2025, total_pop=9_980_999)
sti_pars = dict(hiv=dict(beta_m2f=0.035, eff_condom=0.95, rel_init_prev=1.0))
nw_pars = dict(prop_f0=0.79, prop_m0=0.83, f1_conc=0.16, m1_conc=0.11, p_pair_form=0.58)

def make_sim(**kwargs):
    intvs = make_custom_interventions()
    user_intvs = sc.tolist(kwargs.pop('interventions', []))
    sim = sti.Sim(
        demographics='zimbabwe',
        diseases='hiv',
        data_path=sc.thispath(),
        sim_pars=sim_pars, nw_pars=nw_pars, sti_pars=sti_pars,
        interventions=intvs + user_intvs,
        **kwargs,
    )
    return sim
```

CSV files in the folder:
- `init_prev_hiv.csv` — initial prevalence by risk group/sex/SW status
- `condom_use.csv` — condom use by partnership type over time
- `art_coverage.csv` — ART coverage (n_art by year)
- `vmmc_coverage.csv` — VMMC coverage (n_vmmc by year)
- `zimbabwe_hiv_data.csv` — UNAIDS estimates for plot overlay

When plotted (`sim.plot()`), results show with UNAIDS data points overlaid automatically (via `sim.data`).

### Comparison

| | `simple` | `zimbabwe` |
|---|---|---|
| CSV data files | None | 5 files |
| Demographics | hivsim defaults (`ss.Pregnancy`, `ss.Deaths`) | Real Zimbabwe data from `stisim/data/files/` |
| Disease pars | hivsim defaults | Calibrated values |
| Network pars | Defaults | Calibrated values |
| Interventions | hivsim defaults (basic testing, ART, VMMC, PrEP) | Custom testing scale-up + data-driven ART/VMMC |
| Plot data overlay | No | Yes (UNAIDS estimates) |
| Purpose | Quickstart, show minimal working sim | Realistic country model, template for new projects |

### hivsim.demo() implementation

Added to `hivsim/sim.py`:

```python
def demo(example=None, run=True, plot=True, list=False, **kwargs):
    """
    Create a demo HIVsim simulation.

    Args:
        example (str): Example name ('simple', 'zimbabwe'). Default: 'simple'.
        run (bool): Whether to run the sim.
        plot (bool): Whether to plot results (only if run=True).
        list (bool): If True, print available examples and return.
        **kwargs: Passed to the example's make_sim().

    Returns:
        Sim: Configured (and optionally run) simulation.

    Examples:
        hs.demo()                         # Run simple default demo
        hs.demo('zimbabwe')               # Run Zimbabwe HIV model
        sim = hs.demo('zimbabwe', run=False, n_agents=500)
    """
    import hivsim_examples as hx

    if list:
        hx.list_examples()
        return

    if example is None:
        example = 'simple'

    make_sim = hx.get_example(example)
    sim = make_sim(**kwargs)

    if run:
        sim.run()
        if plot:
            sim.plot()

    return sim
```

### hivsim_examples/__init__.py

```python
import importlib

EXAMPLES = {
    'simple':   'Minimal HIV sim — no data files, hivsim defaults',
    'zimbabwe': 'Zimbabwe HIV model with calibrated parameters and UNAIDS data',
}

def list_examples():
    """Print available examples."""
    for name, desc in EXAMPLES.items():
        print(f'  {name:12s} {desc}')

def get_example(name):
    """Return the make_sim function for a named example."""
    if name not in EXAMPLES:
        available = ', '.join(EXAMPLES.keys())
        raise ValueError(f"Example '{name}' not found. Available: {available}")
    mod = importlib.import_module(f'hivsim_examples.{name}.sim')
    return mod.make_sim
```

### Relationship to external repos

The external project repos (hiv_kenya, hiv_zambia, hiv_zim, stisim_demo_zambia, stisim_vddx_zim) follow the same `make_sim()` + `data/` pattern as the zimbabwe example. Users building their own country model should:

1. Look at `hivsim_examples/zimbabwe/` as the template
2. Look at the external repos for more advanced examples (calibration, partner notification, multi-disease, scenarios)
3. Copy the structure into their own repo

The README in hivsim_examples should point to these external repos as real-world references.

## Implementation Steps

### Step 1: Create hivsim_examples package

- [ ] `hivsim_examples/__init__.py` — EXAMPLES registry, `list_examples()`, `get_example()`
- [ ] `hivsim_examples/simple/__init__.py`
- [ ] `hivsim_examples/simple/sim.py` — `make_sim()` using `hs.Sim()` with no data files
- [ ] `hivsim_examples/zimbabwe/__init__.py`
- [ ] `hivsim_examples/zimbabwe/sim.py` — `make_sim()` with calibrated Zimbabwe pars + custom interventions (from hiv-org branch)
- [ ] `hivsim_examples/zimbabwe/*.csv` — Zimbabwe data files (from hiv-org branch)

### Step 2: Add hivsim.demo()

- [ ] Add `demo()` function to `hivsim/sim.py`
- [ ] Export from `hivsim/__init__.py`

### Step 3: Tests

- [ ] Update `tests/test_examples.py` with tests for the simplified API
- [ ] Test `hs.demo()` and `hs.demo('zimbabwe')` run end-to-end
- [ ] Test `make_sim()` directly with kwarg overrides
- [ ] Test Zimbabwe sim produces plot with data overlay

### Step 4: Documentation

- [ ] `hivsim_examples/README.md` — usage guide, file format specs, "Creating Your Own Country Model" section, links to external repos

## Future: Automated Data Loading (Problem 2)

When `stisim_data` is ready, the integration point is clear:

```python
# This path (currently raises NotImplementedError in loaders.py):
hs.Sim(location='kenya')
# Would auto-download demographic + disease data from stisim_data,
# then use DataLoader exactly as the examples do.
```

The example project structure (CSVs + make_sim) is forward-compatible — the same CSV formats used in examples will be what stisim_data produces.
