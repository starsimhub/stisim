# Calibration

STIsim models are fitted to data with Starsim's `sim.fit()`. The [calibration tutorial](../tutorials/tut_calibration.qmd) walks through a complete example end-to-end; this page summarizes the parts.

## Overview

A calibration needs two things:

1. **Calibration parameters**: what to vary and over what range, given as `ss.par(value, bounds=[low, high])` in place of the parameter's value.
2. **Data**: a dataframe with a `time` or `year` column and a column named after each result to fit, attached to the sim with `data=`.

```python
sim = sti.Sim(
    diseases='hiv',
    hiv=dict(beta_m2f=ss.par(0.05, bounds=[0.02, 0.10])),
    data=data_df,
)
fit = sim.fit(
    weights={'hiv.prevalence': 1.0, 'hiv.on_art': 0.5},
    n_trials=200,
    n_workers=4,
)
```

## Calibration parameters

Any parameter of any module can be fitted: a disease (`sti.HIV(beta_m2f=ss.par(...))`), the network (`sti.StructuredSexual(prop_f0=ss.par(...))`, or `prop_f0=ss.par(...)` given to `sti.Sim`), a connector, or an intervention. Each needs `bounds`; the value given is the starting point. A parameter whose name another fitted module also has (e.g. `beta_m2f` of HIV and of syphilis) is named after its module in the fit (e.g. `syph.beta_m2f`).

A transformation of the parameters (e.g. one calibration parameter that scales several sim parameters) is written as a module that reads it, or by fitting the underlying parameters directly.

## Goodness of fit

By default, the fit compares each result in the data with the corresponding column (a flow summed over each data period, a stock or ratio averaged over it), computes the normalized absolute error, applies `weights`, and sums. Other choices:

| Option | What it does |
|--------|--------------|
| `gof='nae'` | Normalized absolute error (the default) |
| `gof='rmse'` | Root mean squared error |
| `gof=fn` | A function of the sim that returns the loss, e.g. one that returns `np.inf` for implausible runs and `ss.fit.mismatch(sim)` otherwise |
| `likelihood=dict(...)` | A likelihood for each data series (e.g. `ss.binom(n=...)`), for a fit that maximizes the likelihood |

## Running the calibration

`sim.fit()` runs Optuna trials (in parallel with `n_workers=`), and returns a fit with the best values as `fit.pars` and the calibrated sim, already run, as `fit.sim`. After calibration:

```python
fit.summarize()               # Goodness of fit before and after, and the estimates
fit.plot()                    # Each result with data, before and after
sim2 = sim.new(pars=fit.pars) # A new sim with the best-fit values
trials = fit.to_df()          # Every trial: the values and the score (objective)
```

With `storage='fit.log'`, the trials are kept in a file, and calling the fit again resumes it.

## Running the best parameter sets

```python
par_sets = fit.to_df().nsmallest(50, 'objective').drop(columns='objective').to_dict('records')
msim = ss.parallel(fit.sim, pars=par_sets, n_runs=2)
```

Each parameter set is a scenario of its own in the results (`msim.results`), run `n_runs` times with different seeds.

## API reference

- `sim.fit()`: fit the parameters with bounds to the data (see the [Starsim documentation](https://docs.starsim.org))
- `ss.fit.mismatch()`: the default goodness of fit of a sim that has run
- `ss.parallel()`: run sims with sets of parameter values
