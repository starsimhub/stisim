# Data guide

The data an STIsim analysis needs depends on the study type. Route by
that first, then follow the specific path.

## Study type

**Country study** — the analysis targets a real country's HIV / STI
dynamics for policy or projection. Needs demographic data, initial
prevalence, calibration targets, and behavioural inputs. This is the
most common shape.

**Theoretical or illustrative study** — the analysis is about
mechanism or behaviour (transmission patterns, network dynamics,
intervention shapes) rather than a specific setting. Uses default
STIsim demographics or a synthetic population. No country-specific
data files needed. Populations are typically small (5k–20k agents)
and time horizons short (1–5 years).

**Sub-national study** — the analysis targets a subregion (a province,
a district, a key population). Requires location-specific data the
user provides. STIsim's default demographic loaders won't cover it;
the user brings CSVs in the expected format (see below).

## Country study — data requirements

### Demographics

Required to populate the background population dynamics.

| What | Source | Notes |
|---|---|---|
| Age structure by sex (initial year) | UN World Population Prospects (WPP) | `stisim.data.downloaders` fetches this from the UN Data Portal API given an auth key (see below). |
| Age-specific fertility rate (ASFR) | UN WPP or DHS | Available via the same downloader path; DHS is an alternative for a specific recent survey year. |
| Death rates by age and sex | UN WPP life tables | Available via the same downloader path. |
| Net migration | UN WPP or country statistics office | Available via the downloader; may need `rel_migration` scaling to match a specific dataset (example: `dem_pars={'rel_migration': 0.5}`). |

**STIsim ships a UN Data Portal API downloader.** See
`stisim/data/downloaders.py`. To use it:

1. Register for an auth token by emailing `population@un.org` with the
   subject "Data Portal Token Request"
   (`population.un.org/dataportalapi/index.html`).
2. Save the token as `stisim/data/files/auth_key.txt` in your STIsim
   install.
3. Call the downloader for the target location; it fetches death
   rates, ASFR, and (optionally) births, writes them to
   `stisim/data/files/` as `<location>_deaths.csv`,
   `<location>_asfr.csv`, `<location>_births.csv`.

If no auth key is set, the downloader raises an error with manual-
download instructions. Users then download the CSVs by hand from UN
WPP and place them in the same folder with the same filenames.

Expected filenames follow the `<location>_<indicator>.csv` convention
so that `sti.Sim(demographics='<country>')` resolves to the files
without extra configuration. See `tests/test_data/zimbabwe_asfr.csv`,
`zimbabwe_deaths.csv`, `zimbabwe_births.csv` for canonical shape.

**High-HIV-burden countries: deduplicate the deaths file before
running.** UN WPP all-cause mortality includes AIDS deaths. Because
the HIV module separately kills agents via `p_hiv_death` and the
`ti_zero` pathway, feeding raw all-cause rates into `ss.Deaths`
double-counts HIV mortality. Use `stisim.data.dedup_deaths` on the
deaths file before running:

```python
import pandas as pd
from stisim.data import dedup_deaths, deleted_fraction

all_cause = pd.read_csv('data/<country>_deaths.csv')
hiv_deleted = dedup_deaths(all_cause, base_year=1985, end_year=2025)

# Preserve the original under a name stisim will not pick up; write
# the HIV-deleted version to the file stisim reads.
all_cause.to_csv('data/<country>_deaths_all_cause.csv', index=False)
hiv_deleted.to_csv('data/<country>_deaths.csv', index=False)

# Sanity check: AIDS share of all-cause mortality should peak in the
# mid-30s (~70–85%) and fall to ~0 at 80+.
diag = deleted_fraction(all_cause, hiv_deleted)
```

See `model-primer/references/calibration-knobs.md` under "HIV
mortality" for the argument and the full mortality-adjustment
sequence.

### Initial disease prevalence

Required for each disease in the composition.

| Disease | Typical source |
|---|---|
| HIV | National surveillance (Spectrum estimates), PHIA-family surveys (ZIMPHIA, MPHIA, etc.), IBBS for key populations. |
| Syphilis | PHIA-family surveys where included (syphilis is in ZIMPHIA), ANC surveillance, national STI programme data. |
| Gonorrhoea / chlamydia / trichomoniasis | Sparse. Etiologic surveys where they exist; otherwise regional or global estimates (GBD). |

Prevalence files are typically per-disease CSVs (e.g. `init_prev_hiv.csv`)
with age × sex breakdowns matching STIsim's expected age bins (see
`stisim/diseases/sti.py` `default_age_bins`).

### Calibration targets

Time series the analysis will calibrate against.

| Source | What it provides |
|---|---|
| IHME Global Burden of Disease | Country-year prevalence / incidence / mortality for HIV and the STIs. Raw files are large and typically not committed to the repo. See `stisim_vddx_zim/process_ihme_data.py` for a local-processing pattern. |
| UNAIDS Spectrum | HIV-specific: prevalence, ART coverage, incidence, mortality. |
| National surveillance | Any country-specific programme data — new HIV diagnoses, syphilis notifications, ANC positivity. |

### Behavioural inputs

| Input | Source |
|---|---|
| Condom use by partnership type and year | DHS (Domestic Health Survey). Typical file: `condom_use.csv`. |
| Age at sexual debut, partnership formation rates | DHS, if the analysis needs to move them off defaults. |
| Sex-worker population parameters | IBBS surveys, national key population size estimates. |

## Sub-national study — data expectations

The user provides all location-specific data. STIsim's default loaders
won't fetch subregional data. Files should follow the same shape as
the country CSVs above:

- Demographic CSVs with columns `age`, `sex`, `year`, `value` (units
  as appropriate: population count for age structure, per-1000-women
  for ASFR, per-1000-population for mortality).
- Prevalence CSVs with `age`, `sex`, `value` (proportion).
- Calibration target CSVs with `year`, `value`, `lower`, `upper`
  (bounds for CI-based likelihood).

Point STIsim at the user's `data/` directory rather than using the
`demographics='<country>'` string convention — the string convention
resolves against STIsim's packaged data loaders, not user files.

## Theoretical study — no external data

Populations of 5k–20k agents run over 1–5 years illustrate transmission
dynamics without needing country data. Skip data ingestion entirely.
Use `sti.Sim()` with default demographics and small `n_agents`. Fix
seeds and vary only the axis of interest (network structure,
intervention shape) across runs.

## File layout

Convention across the surveyed projects:

```
<project-root>/
├── data/
│   ├── <country>_age_male.csv
│   ├── <country>_age_female.csv
│   ├── <country>_asfr.csv
│   ├── <country>_deaths.csv
│   ├── init_prev_<disease>.csv
│   ├── <country>_sti_data.csv          # calibration targets
│   ├── <country>_hiv_data.csv          # calibration targets
│   └── condom_use.csv
├── data-raw/                           # optional; not committed
│   └── ...                             # large IHME extracts etc.
└── ...
```

Files derived from large raw sources (IHME) typically go through a
`process_*.py` script that's committed alongside the derived CSV; the
raw extract itself is left out of the repo.
