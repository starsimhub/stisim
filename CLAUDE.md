# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

STIsim is an agent-based modeling framework for sexually-transmitted infections, built on top of Starsim (starsim). It provides disease models (HIV, syphilis, chlamydia, gonorrhea, trichomonas, BV, GUD), structured sexual networks with risk groups, demographics, and interventions.

## Build & Test Commands

```bash
# Install in development mode
pip install -e .

# Run all tests (from tests/ directory)
cd tests && pytest test_*.py -n auto

# Run a specific test
pytest tests/test_sim.py::test_hiv_sim -v

# Run tests matching a pattern
pytest tests/test_*.py -k "hiv" -v

# Run with coverage
pytest tests/test_*.py --cov=stisim
```

Tests use pytest with xdist for parallel execution. The pytest.ini sets `SCIRIS_BACKEND=agg` for non-interactive plotting.

## Architecture

**Inheritance**: STIsim extends Starsim (`ss.*` classes). The main `Sim` class inherits from `ss.Sim`.

**Module hierarchy**:
- `stisim/sim.py` - Core Sim class, orchestrates simulations
- `stisim/diseases/` - Disease modules (each extends base classes in `sti.py`)
- `stisim/networks.py` - Sexual contact networks (StructuredSexual, PriorPartners, AgeMatchedMSM)
- `stisim/interventions/` - Testing, treatment, partner notification
- `stisim/connectors/` - Disease coinfection interactions (HIV-STI, GUD-syphilis)
- `stisim/analyzers.py` - Result tracking (coinfection_stats, sw_stats, etc.)
- `stisim/demographics.py` - Pregnancy, migration modules
- `stisim/calibration.py` - Model calibration tools

**Key pattern**: Each module has a `*Pars` class for parameters (e.g., `CTPars`, `NetworkPars`, `SimPars`).

## Code Conventions

**Imports**:
```python
import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
```

**Naming**:
- Parameter classes: suffix with `Pars` (e.g., `HIVPars`)
- State booleans: prefix with `is_` or use action names (tested, diagnosed, cleared)
- Time intervals: prefix with `ti_` (e.g., `ti_infected`, `ti_diagnosed`)

**Style**: Follows Google Python style guide (https://github.com/starsimhub/styleguide). No automated linting in CI.

## Common Patterns

**Creating a simulation**:
```python
sim = sti.Sim(
    diseases=['hiv', 'syphilis'],  # or disease instances
    networks=sti.StructuredSexual(),
    n_agents=1000,
    dur=40
)
sim.run()
```

**Disease implementation** (in `stisim/diseases/`):
- Extend `SEIS` for standard STIs, or specialized base classes for HIV/syphilis/BV
- Define `__init__`, states, transmission logic
- Use `define_pars()` and `update_pars()` pattern

**State management**:
```python
self.define_states(
    ss.BoolState('infected'),
    ss.FloatArr('ti_infected'),
)
```

**Result tracking**:
```python
self.define_results(
    ss.Result('incidence', dtype=int, scale=True),
    ss.Result('prevalence', dtype=float, scale=False),
)
```

## Key Parameters

**Transmission**: `beta_m2f`, `beta_f2m`, `beta_m2m`, `beta_m2c` (male-to-commercial)

**Natural history**: `dur_exp` (latent), `p_symp` (symptomatic probability), `dur_presymp`, `dur_asymp2clear`

**Networks**: Risk groups (0=low, 1=mid, 2=high), age preferences, concurrency, relationship duration

## Dependencies

Core: `starsim>=3.0.1`, `sciris>=3.1.6`, `pandas`, `scipy`, `numba`, `networkx`, `optuna`

## Resources

- Docs: https://docs.idmod.org/projects/stisim
- Starsim docs: https://docs.starsim.org
- GitHub: https://github.com/starsimhub/stisim
