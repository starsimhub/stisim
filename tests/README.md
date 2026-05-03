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

## Adding tests

- Name new test files `test_*.py` so they are picked up automatically.
- Prefer end-to-end tests that exercise full sim runs.
- When changing baseline behavior, regenerate `baseline.yaml` deliberately and call out the change in the PR description.
