# Data loaders and downloaders

Utilities for fetching and loading demographic and epidemiological inputs used by STIsim simulations.

| File | Contents |
|------|----------|
| `loaders.py` | Functions that load locally cached data (age structures, fertility, mortality, prevalence) into formats consumed by `sti.Sim` and the disease modules. Imported as `stidata` inside `sim.py`. |
| `downloaders.py` | Functions that fetch data from upstream sources (UN World Population Prospects, UNAIDS, etc.) and cache them on disk. Imported as `stidl` inside `sim.py`. |
| `test_downloaders.py` | Smoke tests for the downloaders (kept here rather than in `tests/` because they exercise package-internal paths and are skipped by default in CI). |

Most users do not call these directly — `sti.Sim` invokes the loaders to populate defaults when no data are passed in. Use them explicitly when you need to refresh cached data or build a country-specific scenario from scratch.
