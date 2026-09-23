# AGENTS.md

This file provides guidance to Claude Code (claude.ai/code) and other agents when working with code in this repository.

## Overview

STIsim is an agent-based modeling framework for co-circulating sexually transmitted infections, built on [Starsim](https://docs.starsim.org) (`starsim>=3.6.0`). Nearly every class subclasses a Starsim base (`ss.Sim`, `ss.Infection`, `ss.SexualNetwork`, `ss.Intervention`, `ss.Connector`, `ss.Pars`), so Starsim semantics (module lifecycle, `ss.Arr` states, UIDs, time/rate types like `ss.years`/`ss.permonth`, results) apply throughout. Convention: `import stisim as sti`, `import starsim as ss`, `import sciris as sc`.

The repo ships three installable packages: `stisim` (the library), `hivsim` (a thin HIV-defaults wrapper), and `hivsim_examples` (runnable example sims plus Zimbabwe CSV data).

## Commands

```sh
pip install -e .[dev]                   # Dev install (includes pytest, xdist, quarto/jupyter tooling)

# Tests — run from tests/
cd tests
./run_tests                             # Full suite in parallel; sets SCIRIS_BACKEND=agg for non-interactive plots
pytest test_*.py -n auto                # What CI runs (Python 3.13)
pytest test_sim.py::test_minimal_hiv -v # Single test
pytest test_*.py -k hiv -v              # Pattern match
./update_baseline                       # Regenerate baseline.yaml and benchmark.yaml
```

- `tests/test_baselines.py` runs `hivsim.demo('zimbabwe')` (2000 agents, to 2010, seed 2) and compares its summary against `tests/baseline.yaml` via `ss.diff_sims(..., die=True)`. Any change that alters model dynamics or random-number draws will fail this; regenerate the baseline deliberately with `./update_baseline` and call it out in the PR.
- `tests/devtests/` holds exploratory scripts, not part of the suite. `stisim/data/test_downloaders.py` needs network access and is not collected by CI.
- Docs are Quarto + quartodoc (`docs/_quarto.yml`); from `docs/`: `./render` (full build), `./preview` (live reload). `docs/api/` is generated. Tutorials and examples are `.qmd` files.

## Architecture

### `sti.Sim` construction and parameter routing (`stisim/sim.py`, `stisim/parameters.py`, `stisim/utils.py`)

`sti.Sim` is where most of the non-obvious logic lives. It accepts modules in several forms (strings like `diseases=['hiv', 'syph']`, module instances, or flat/nested pars) and resolves them in two phases:

1. **`__init__` → `separate_pars`**: calls the parent `ss.Sim.__init__` *without* modules, then sorts all incoming pars into sim / sti / network / demographic / connector buckets via `sti.route_pars`. Routing is registry-based: a flat key is matched against the par classes of every registered disease (`sti_register()`), `NetworkPars`, demographic pars (`dem_pars()`), and auto-connectors (`merged_connector_pars()`). Keys shared by several diseases (e.g. `beta_m2f`, `init_prev`) are intentionally broadcast to all of them; per-disease dicts like `hiv=dict(...)` scope pars to one disease. `remap_pars` handles legacy aliases (`beta` → `beta_m2f`, `location` → `demographics`).
2. **`init`**: builds modules from the stashed pars: `process_stis` (string names resolved through `sti_aliases()` / `make_sti`), `process_connectors`, `process_networks` (default: `StructuredSexual` + `ss.MaternalNet`), and `process_demographics`, then calls `ss.Sim.init`.

**XOR rule:** passing a module instance *and* pars targeting that same slot (e.g. `diseases=[sti.HIV()]` plus `init_prev=...`, or `networks=[...]` plus network pars) raises a `ValueError` rather than guessing. The `_user_*_pars` attributes exist to enforce this.

**Adding a new disease or connector** means registering it: diseases in `sti_register()` and `sti_aliases()`; auto-discoverable connectors in `connector_register()`. Otherwise its pars won't route and string construction won't work.

**Location-based demographics:** `demographics='zimbabwe'` (or `location=`) makes `process_demographics` check for cached UN WPP data under `stisim/data`, download it if missing (requires internet), and build `ss.Pregnancy`/`ss.Births`, `ss.Deaths`, optional `sti.Migration`, and `ss.People` from the age distribution. `total_pop`/`age_scale` control scaling.

### Diseases (`stisim/diseases/`)

`BaseSTI(ss.Infection)` in `sti.py` is the common base; `SEIS(BaseSTI)` is the template for bacterial STIs (chlamydia, gonorrhea, trichomoniasis). HIV, syphilis and BV subclass `BaseSTI` directly with their own natural history. Each disease has a companion `*Pars` class (`BaseSTIPars` → `STIPars` → disease-specific). Transmission is expressed per network through `BaseSTI.validate_beta`, which maps `beta_m2f` (+ `rel_beta_f2m`) → `structuredsexual`, `beta_m2c` → `maternal`, `beta_breastfeed` → `breastfeeding`, `beta_m2m` → `msm`. Network *names* therefore matter: a custom network needs one of those names for these pars to apply. Results are age/sex stratified by default; `age_bins=None` / `sex_keys=None` disable stratification.

### Networks (`stisim/networks/`)

`BaseNetwork(ss.SexualNetwork)` → `MFNetwork` → `StructuredSexual` (the default: risk groups, partnership types, sex work, debut ages). Partner matching algorithms live in `matchers.py` and are pluggable. MSM networks are in `msm.py`.

### Interventions (`stisim/interventions/`)

Products (diagnostics with sensitivity/specificity, e.g. `STIDx`) are separate from interventions (eligibility, timing, coverage: `STITest`, `SymptomaticTesting`, `STITreatment`, `PartnerNotification`). HIV interventions (`HIVTest`, `ART`, `VMMC`, `Prep`) are in `hiv_interventions.py`. `stisim/logistics/` adds a supply-chain layer (products, supplies, `SuppliedIntervention`). `care_seeking.py` provides a per-agent cross-disease care-seeking propensity that interventions can use.

### Connectors (`stisim/connectors/`)

Coinfection effects that adjust `rel_sus`/`rel_trans` each step. `process_connectors` auto-adds any `sti.<d1>_<d2>` class (e.g. `hiv_syph`) for each disease pair in the sim, unless the user supplied `connectors=[...]`. Flat connector pars passed to `sti.Sim` are routed only to these auto-added defaults.

### HIVsim (`hivsim/sim.py`)

`hivsim.Sim` is a thin `sti.Sim` subclass that sets HIV defaults (HIV disease; Pregnancy + Deaths; StructuredSexual + Maternal + Breastfeeding networks; HIVTest/ART/VMMC/PrEP) and hands all par routing to `sti.Sim`. `hivsim.demo(name)` loads examples from `hivsim_examples/` (`simple`, `zimbabwe`).

### `stisim.ai` (`stisim/ai/`)

A Claude Code plugin shipped inside the package (skills under `stisim/ai/plugin/skills/`, registered with `python -m stisim.ai install`). Skill Markdown is package data (see `[tool.setuptools.package-data]` in `pyproject.toml`), so new file types under `plugin/` must be added there to ship.

## Conventions

- Version and date live in `stisim/version.py`; the minimum Starsim version is duplicated in `pyproject.toml` and the `sc.require` call in `stisim/__init__.py`, so keep them in sync. User-facing changes go in `CHANGELOG.md` (`docs/whatsnew.md` is a symlink to it).
- Code style follows the [Starsim style guide](https://github.com/starsimhub/styleguide) (roughly Google Python style with exceptions); docstrings are Google-style, since quartodoc parses them for the API reference.
- PRs target `main`.
