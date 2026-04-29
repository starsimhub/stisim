# Contributing

Welcome! We are thrilled you are interested in contributing to STIsim. This document will help you get started.

## Code of conduct and style

- The STIsim community follows a [code of conduct](https://docs.idmod.org/projects/starsim/en/stable/conduct.html). By participating in this project, you agree to abide by its terms.
- Code should follow our house [style guide](https://github.com/starsimhub/styleguide). STIsim more or less follows Google's Python style guide, with some exceptions.

## Development install

Clone the repository and install in editable mode with the development extras:

```sh
git clone https://github.com/starsimhub/stisim.git
cd stisim
pip install -e .[dev]
```

The `[dev]` extras include `pytest`, `pytest-xdist`, and the documentation toolchain (`mkdocs`, `mkdocstrings`, `mkdocs-jupyter`).

## Running the tests

Tests use `pytest` with `pytest-xdist` for parallel execution. From the `tests/` directory:

```sh
cd tests
pytest test_*.py -n auto                 # Run the full test suite in parallel
pytest test_sim.py::test_hiv_sim -v      # Run a specific test
pytest test_*.py -k "hiv" -v             # Run tests matching a pattern
pytest test_*.py --cov=stisim            # Run with coverage
```

`pytest.ini` sets `SCIRIS_BACKEND=agg` for non-interactive plotting. Tests also run automatically on every push and pull request via GitHub Actions.

## Building the documentation

```sh
mkdocs serve                             # Live-reload preview at http://127.0.0.1:8000
mkdocs build                             # Static build into site/
```

The build executes all tutorial notebooks. If you change a public API, run `mkdocs build` locally before opening a pull request to confirm tutorials still execute cleanly.

## Pull requests

- Branch off and target `main`.
- Keep PRs focused — one logical change per PR.
- Add or update tests for any behavior change.
- Add an entry to `CHANGELOG.md` under the next unreleased version, with the PR number.
- Update relevant documentation (user guide, tutorials, docstrings) alongside code changes.
- Make sure `pytest test_*.py -n auto` passes locally before requesting review.

## Issues

Feel free to [open an issue](https://github.com/starsimhub/stisim/issues/new/choose) on more or less anything — bug reports, feature requests, or questions. This project is small enough that we don't need a formal triage system.

## Releasing

See `.github/workflows/README.md` for details on cutting and publishing releases.

## Questions

If you have any other questions, please reach out to us: [info@starsim.org](mailto:info@starsim.org). Thank you!
