# Calibration

STIsim ships a thin wrapper around [Optuna](https://optuna.org/) for fitting model parameters to data. The [calibration tutorial](../tutorials/tut_calibration.ipynb) walks through a complete example end-to-end; this page documents the public surface.

## Overview

A calibration run requires three things:

1. **Calibration parameters** — what to vary and over what range.
2. **A build function** — how to apply a candidate parameter set to a sim.
3. **An evaluation function** — how to score a candidate against data.

```python
calib = sti.Calibration(
    sim=base_sim,
    calib_pars={'hiv.beta_m2f': dict(low=0.02, high=0.10, guess=0.05)},
    data=data_df,
    sim_result_list=['hiv.prevalence', 'hiv.n_on_art'],
    weights={'hiv.prevalence': 1.0, 'hiv.n_on_art': 0.5},
    total_trials=200,
    n_workers=4,
)
calib.calibrate()
```

## Calibration parameters

`calib_pars` is a flat dict keyed by `'<module>.<parameter>'`. Each entry is itself a dict with `low`, `high`, and `guess` keys, optionally `suggest_type` (`'float'`, `'int'`, `'categorical'`).

Keys are resolved via `flatten_calib_pars()` and applied via `set_sim_pars()` — these are exposed for use in custom build functions.

## Build function

The default build function applies parameters in place. Override it when calibration parameters need to be transformed before being applied (e.g., a single calibration parameter that scales several sim parameters):

```python
def my_build_fn(sim, calib_pars, **kwargs):
    factor = calib_pars.pop('beta_scale')
    sim = sti.set_sim_pars(sim, calib_pars)
    for disease in sim.diseases.values():
        if hasattr(disease.pars, 'beta_m2f'):
            disease.pars.beta_m2f *= factor
    return sim
```

## Evaluation function

The default evaluation function computes a goodness-of-fit metric for each result in `sim_result_list` against the corresponding column in `data`, applies `weights`, and returns the weighted sum.

`compute_gof()` supports several metrics via the `as_scalar` argument:

| `as_scalar` | Metric |
|-------------|--------|
| `'none'` | Element-wise normalized errors (returns array) |
| `'sum'` | Sum of normalized errors |
| `'mean'` | Mean normalized error |
| `'median'` | Median normalized error (robust) |

For custom likelihoods, supply your own `eval_fn(sim, data=...)` callable.

## Running the calibration

`calib.calibrate()` runs Optuna trials in parallel and stores the best parameters on `calib.best_pars`. After calibration:

```python
sim = sti.set_sim_pars(base_sim, calib.best_pars)
sim.run()
```

`calib.plot()` provides parameter-importance and trial-history plots.

## API reference

- `sti.Calibration` — main class
- `sti.compute_gof` — goodness-of-fit metrics
- `sti.flatten_calib_pars` — flatten nested calib pars to dotted keys
- `sti.set_sim_pars` — apply a flat dict of parameters to a sim
- `sti.default_build_fn` — default build function
- `sti.eval_fn` — default evaluation function
