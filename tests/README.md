# Tests

STIsim's test suite. Tests use `pytest` with `pytest-xdist` for parallel execution.

## Running

From this directory:

```sh
pytest test_*.py -n auto                 # Full suite, parallel
pytest test_sim.py::test_hiv_sim -v      # Single test
pytest test_*.py -k "hiv" -v             # Tests matching a pattern
pytest test_*.py --cov=stisim            # With coverage
./run_tests                              # Convenience wrapper
```

`pytest.ini` sets `SCIRIS_BACKEND=agg` so plotting is non-interactive. The same suite runs on every push and pull request via GitHub Actions.

## Layout

| File / folder | Purpose |
|---|---|
| `test_sim.py` | Top-level `sti.Sim` construction, parameter routing, and run sanity checks. |
| `test_hiv.py` | HIV disease dynamics (CD4 progression, transmission, MTCT). |
| `test_hiv_interventions.py` | HIV intervention behavior (testing, ART, VMMC, PrEP). |
| `test_stis.py` | Bacterial/viral STI modules (chlamydia, gonorrhea, trichomoniasis, syphilis, BV, GUD). |
| `test_networks.py` | Network construction and partnership formation. |
| `test_connectors.py` | Coinfection connector logic. |
| `test_calibration.py` | Calibration helpers and end-to-end fitting. |
| `test_baselines.py` | Regression tests against committed baseline output. |
| `baseline.yaml`, `benchmark.yaml` | Committed baseline results and runtime benchmarks consumed by `test_baselines.py`. |
| `hiv_natural_history_analyzers.py` | Custom analyzers used by HIV tests; not a test module itself. |
| `simple.py` | Minimal example referenced from tests as a smoke-test target. |
| `testlib.py` | Shared test fixtures and helpers. |
| `test_data/` | Test fixtures (CSVs and small inputs). |
| `devtests/` | Developer-only experiments and exploratory scripts; not run by CI. |

## Adding tests

- Name new test files `test_*.py` so they are picked up automatically.
- Prefer small, fast unit-style tests; reserve longer end-to-end checks for `test_baselines.py`.
- When changing baseline behavior, regenerate `baseline.yaml` deliberately and call out the change in the PR description.
