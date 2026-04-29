"""
Define the calibration class
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import os

# Lazy import (do not import unless actually used, saves load time)
op = sc.importbyname('optuna', lazy=True)

__all__ = ['Calibration', 'compute_gof', 'default_build_fn', 'set_sim_pars', 'flatten_calib_pars', 'make_calib_sims']


def compute_gof(actual, predicted, normalize=True, use_frac=False, use_squared=False,
                as_scalar='none', eps=1e-9, skestimator=None, estimator=None, **kwargs):
    """
    Calculate the goodness of fit. By default use normalized absolute error, but
    highly customizable. For example, mean squared error is equivalent to
    setting normalize=False, use_squared=True, as_scalar='mean'.

    Args:
        actual      (arr):   array of actual (data) points
        predicted   (arr):   corresponding array of predicted (model) points
        normalize   (bool):  whether to divide the values by the largest value in either series
        use_frac    (bool):  convert to fractional mismatches rather than absolute
        use_squared (bool):  square the mismatches
        as_scalar   (str):   return as a scalar instead of a time series: choices are sum, mean, median
        eps         (float): to avoid divide-by-zero
        skestimator (str):   if provided, use this scikit-learn estimator instead
        estimator   (func):  if provided, use this custom estimator instead
        kwargs      (dict):  passed to the scikit-learn or custom estimator

    Returns:
        gofs (arr): array of goodness-of-fit values, or a single value if as_scalar is True

    **Examples**::

        x1 = np.cumsum(np.random.random(100))
        x2 = np.cumsum(np.random.random(100))

        e1 = compute_gof(x1, x2) # Default, normalized absolute error
        e2 = compute_gof(x1, x2, normalize=False, use_frac=False) # Fractional error
        e3 = compute_gof(x1, x2, normalize=False, use_squared=True, as_scalar='mean') # Mean squared error
        e4 = compute_gof(x1, x2, skestimator='mean_squared_error') # Scikit-learn's MSE method
        e5 = compute_gof(x1, x2, as_scalar='median') # Normalized median absolute error -- highly robust
    """

    # Handle inputs
    actual    = np.array(sc.dcp(actual), dtype=float)
    predicted = np.array(sc.dcp(predicted), dtype=float)

    # Scikit-learn estimator is supplied: use that
    if skestimator is not None: # pragma: no cover
        try:
            import sklearn.metrics as sm
            sklearn_gof = getattr(sm, skestimator) # Shortcut to e.g. sklearn.metrics.max_error
        except ImportError as E:
            errormsg = f'You must have scikit-learn >=0.22.2 installed: {str(E)}'
            raise ImportError(errormsg) from E
        except AttributeError as E:
            errormsg = f'Estimator {skestimator} is not available; see https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter for options'
            raise AttributeError(errormsg) from E
        gof = sklearn_gof(actual, predicted, **kwargs)
        return gof

    # Custom estimator is supplied: use that
    if estimator is not None: # pragma: no cover
        try:
            gof = estimator(actual, predicted, **kwargs)
        except Exception as E:
            errormsg = f'Custom estimator "{estimator}" must be a callable function that accepts actual and predicted arrays, plus optional kwargs'
            raise RuntimeError(errormsg) from E
        return gof

    # Default case: calculate it manually
    else:
        # Key step -- calculate the mismatch!
        gofs = abs(np.array(actual) - np.array(predicted))

        if normalize and not use_frac:
            actual_max = abs(actual).max()
            if actual_max > 0:
                gofs /= actual_max

        if use_frac:
            if (actual<0).any() or (predicted<0).any():
                print('Warning: Calculating fractional errors for non-positive quantities is ill-advised!')
            else:
                maxvals = np.maximum(actual, predicted) + eps
                gofs /= maxvals

        if use_squared:
            gofs = gofs**2

        if as_scalar == 'sum':
            gofs = np.sum(gofs)
        elif as_scalar == 'mean':
            gofs = np.mean(gofs)
        elif as_scalar == 'median':
            gofs = np.median(gofs)

        return gofs


def make_df(sim, df_res_list=None):
    """
    Build a year-resampled DataFrame of selected results from a sim or MultiSim.

    Used by `eval_fn` to align modeled output with observed data on a
    common yearly time axis before computing goodness of fit.

    Args:
        sim (Sim or MultiSim): Completed simulation(s) to extract results from.
        df_res_list (list[str]): Result keys to include, either at the top level
            (`'hiv.prevalence'`) or as `'<module>.<result>'` paths into a
            nested `Results` container.

    Returns:
        pandas.DataFrame: One column per requested result plus a `'time'`
        column (year). Rows are stacked across sims when `sim` is a MultiSim.
    """
    dfs = sc.autolist()

    # sometimes sim is a MultiSim, so we need to handle that
    if isinstance(sim, ss.MultiSim):
        df_sims = sim.sims
    else:
        df_sims = [sim]

    for df_sim in df_sims:
        for sres in df_res_list:
            if sres in df_sim.results.keys() and isinstance(df_sim.results[sres], ss.Result):
                dfs += df_sim.results[sres].to_df(resample='year', use_years=True, col_names=sres)
            else:
                rescomps = sres.split('.')
                modname = rescomps[0]
                if isinstance(df_sim.results[modname], ss.Results):
                    resname = '_'.join(rescomps[1:])
                    dfs += df_sim.results[modname][resname].to_df(resample='year', use_years=True, col_names=sres)

    df_res = pd.concat(dfs, axis=1)
    # Remove duplicate timevec columns from concat, then extract year
    if 'timevec' in df_res.columns:
        timevec = df_res['timevec'].iloc[:, 0] if isinstance(df_res['timevec'], pd.DataFrame) else df_res['timevec']
        df_res = df_res.loc[:, ~df_res.columns.duplicated()]
        df_res['time'] = timevec.dt.year
    else:
        df_res['time'] = df_res.index
    return df_res


def eval_fn(sim, data=None, sim_result_list=None, weights=None, df_res_list=None):
    """
    Default evaluation function for STIsim calibration.

    For each requested result, merges modeled output with observed data on
    `time` (year), computes per-point goodness-of-fit via
    `compute_gof`, applies an optional weight, and sums the result.

    Args:
        sim (Sim): Completed simulation. The intermediate DataFrame is also
            stored as `sim.df_res` for inspection.
        data (dict): Mapping of result key → DataFrame indexed by year, with
            an observed-value column matching the result key.
        sim_result_list (list[str]): Result keys to evaluate (e.g.
            `['hiv.prevalence', 'hiv.n_on_art']`).
        weights (dict): Optional per-result weight; defaults to 1.0.
        df_res_list (list[str]): Result keys to extract via `make_df`;
            defaults to `sim_result_list`.

    Returns:
        float: Total weighted mismatch (lower is better).
    """
    df_res = make_df(sim, df_res_list=df_res_list)
    sim.df_res = df_res

    # Compute fit
    fit = 0
    for skey in sim_result_list:
        data_df = data[skey].reset_index()
        model_df = df_res[['time', skey]]

        combined = pd.merge(data_df, model_df, how='left', on='time')
        combined['diffs'] = combined[skey+'_x'] - combined[skey+'_y']
        gofs = compute_gof(combined.dropna()[skey+'_x'], combined.dropna()[skey+'_y'])
        losses = gofs
        if weights and (skey in weights.keys()) and (weights[skey] != 1):
            losses *= weights[skey]
        mismatch = losses.sum()
        fit += mismatch

    return fit


# Metadata keys to skip when routing parameters
_META_KEYS = {'rand_seed', 'index', 'mismatch'}
# Keys that indicate a parameter spec dict (not a nested module dict)
_SPEC_KEYS = {'low', 'high', 'guess', 'value', 'suggest_type'}


def flatten_calib_pars(calib_pars):
    """
    Normalize calibration parameters to flat dot-notation format.

    Accepts two formats:

    **Nested** (grouped by module)::

        dict(hiv=dict(beta_m2f=dict(low=0.01, high=0.10)))

    **Flat** (dot notation — returned unchanged)::

        {'hiv.beta_m2f': dict(low=0.01, high=0.10)}

    Nested format is detected when a value is a dict whose keys don't
    overlap with spec keys (``low``, ``high``, ``guess``, etc.).

    Returns:
        dict: Flat dict with ``'module.par'`` keys.
    """
    flat = {}
    for key, value in calib_pars.items():
        if not isinstance(value, dict) or (set(value.keys()) & _SPEC_KEYS):
            flat[key] = value  # Already flat, scalar, or metadata
        else:
            for par_name, spec in value.items():
                flat[f'{key}.{par_name}'] = spec
    return flat


def _parse_calib_key(key):
    """
    Parse a calibration parameter key into (module_name, par_name).

    Uses dot notation: ``'hiv.beta_m2f'`` → ``('hiv', 'beta_m2f')``.
    Keys without a dot or in :data:`_META_KEYS` return ``(None, None)``.
    """
    if key in _META_KEYS:
        return None, None
    if '.' not in key:
        return None, None
    mod_name, par_name = key.split('.', 1)
    return mod_name, par_name


def _parse_legacy_key(key, sim):
    """
    Parse a legacy underscore-separated key (e.g. ``'hiv_beta_m2f'``) by
    greedy-matching against the sim's known module names (longest first).

    Also supports common aliases: ``nw_`` → ``structuredsexual``,
    ``conn_`` → first connector, ``syph_`` → ``syph``.

    Returns ``(module_name, par_name)`` or ``(None, None)`` if unmatched.
    """
    # Collect all module names, plus common aliases
    mod_names = [m.name for m in sim.get_modules()]
    aliases = {}
    for name in mod_names:
        if 'structuredsexual' in name:
            aliases['nw'] = name
    for m in sim.get_modules():
        if isinstance(m, ss.Connector):
            aliases['conn'] = m.name

    all_names = list(aliases.items()) + [(n, n) for n in mod_names]
    # Sort by prefix length descending (greedy match)
    all_names.sort(key=lambda x: len(x[0]), reverse=True)

    for prefix, mod_name in all_names:
        if key.startswith(prefix + '_'):
            par_name = key[len(prefix) + 1:]
            if par_name:
                return mod_name, par_name
    return None, None


def _set_par(module, par_name, value):
    """Set a parameter on a module, calling ``.set()`` for distributions."""
    par = module.pars[par_name]
    if hasattr(par, 'set'):
        par.set(value)
    else:
        module.pars[par_name] = value


def set_sim_pars(sim, pars):
    """
    Set calibrated parameters on a sim.

    All parameters are set directly on the module objects via
    ``sim.get_module()``.  This works on uninitialized sims because every
    module stores its ``pars`` dict immediately on construction.

    Supports both dot notation (``'hiv.beta_m2f'``) and legacy underscore
    format (``'hiv_beta_m2f'``).  Legacy keys are matched greedily against
    the sim's module names (longest match first).

    Args:
        sim (Sim): A simulation (modules must be instances, not strings)
        pars (dict): Flat parameter dict, e.g. ``{'hiv.beta_m2f': 0.05, ...}``

    Returns:
        Sim: The same sim, modified in place
    """
    for key, value in pars.items():
        mod_name, par_name = _parse_calib_key(key)

        # Fall back to legacy underscore matching if no dot found
        if mod_name is None and key not in _META_KEYS:
            mod_name, par_name = _parse_legacy_key(key, sim)

        if mod_name is None:
            continue
        if isinstance(value, dict):
            value = value.get('value', None)
        if value is None:
            continue
        module = sim.get_module(mod_name, die=False)
        if module is None:
            print(f'Warning: module "{mod_name}" not found for parameter "{key}"')
            continue
        _set_par(module, par_name, value)
    return sim


def default_build_fn(sim, calib_pars, **kwargs):
    """
    Default build function for STIsim calibration.

    Routes calibration parameters to the correct sim module using dot notation:

        - ``'hiv.beta_m2f'``              → ``sim.get_module('hiv').pars['beta_m2f']``
        - ``'structuredsexual.prop_f0'``   → ``sim.get_module('structuredsexual').pars['prop_f0']``
        - ``'hiv_syph.rel_sus_syph_hiv'``  → ``sim.get_module('hiv_syph').pars['rel_sus_syph_hiv']``
        - ``'symp_algo.rel_test'``         → ``sim.get_module('symp_algo').pars['rel_test']``

    All parameters are set on the uninitialized sim before ``sim.init()`` is
    called.  This works because every module stores its ``pars`` dict
    immediately on construction.

    ``rand_seed`` is handled specially: ``set_sim_pars`` skips it (it is in
    ``_META_KEYS``), so it is applied directly to ``sim.pars`` here.  This
    ensures that ``reseed=True`` in :class:`Calibration` actually uses a
    different seed per trial, and that the stored ``rand_seed`` column in
    ``calib.df`` can be used to reproduce each trial exactly.

    Args:
        sim (Sim): An uninitialized simulation (modules must be instances, not strings)
        calib_pars (dict): Calibration parameters with values set by the sampler

    Returns:
        Sim: The initialized and modified simulation

    Example::

        calib_pars = {
            'hiv.beta_m2f':               dict(low=0.01, high=0.10, guess=0.035, value=0.05),
            'structuredsexual.f1_conc':    dict(low=0.005, high=0.3, guess=0.16, value=0.10),
            'hiv_syph.rel_sus_syph_hiv':   dict(low=1.0, high=4.0, guess=2.5, value=3.0),
            'symp_algo.rel_test':          dict(low=0.5, high=1.5, guess=1.0, value=1.2),
        }
        sim = default_build_fn(sim, calib_pars)
    """
    set_sim_pars(sim, calib_pars)

    # Apply rand_seed explicitly — set_sim_pars skips it (in _META_KEYS),
    # but reseed=True expects each trial to run with its own seed.
    if calib_pars is not None and 'rand_seed' in calib_pars:
        sim.pars['rand_seed'] = int(calib_pars['rand_seed'])

    if not sim.initialized:
        sim.init()

    return sim


def make_calib_sims(calib=None, calib_pars=None, sim=None, n_parsets=None,
                    check_fn=None, seeds_per_par=1, **kwargs):
    """
    Create and run simulations using calibrated parameters.

    Provide either a :class:`Calibration` object or a set of calibration
    parameters directly. The function creates one sim per parameter set (with
    optional seed replication), runs them in parallel, and optionally filters
    the results.

    Args:
        calib (Calibration): A completed calibration. Extracts pars via
            ``calib.get_pars()`` and uses ``calib.sim`` as the base if *sim*
            is not provided.
        calib_pars: Parameter source — one of:

            - **dict**: single parameter set → 1 sim (× *seeds_per_par*)
            - **list[dict]**: N parameter sets → N sims
            - **DataFrame**: rows are parameter sets (like ``calib.df``)

        sim (Sim): Base (uninitialized) simulation to copy. If ``None``, uses
            ``calib.sim``.
        n_parsets (int): Number of top parameter sets to use. ``None`` = all.
        check_fn (callable): Post-run filter — ``check_fn(sim) → bool``.
            Sims returning ``False`` are dropped. If ``None`` and *calib* is
            provided, uses ``calib.check_fn``.
        seeds_per_par (int): Random seeds per parameter set. When > 1, each
            par set is run with multiple seeds and only the first surviving
            seed (per ``check_fn``) is kept.
        **kwargs: Passed to ``ss.parallel()``.

    Returns:
        ss.MultiSim: A MultiSim containing the completed simulations.

    Examples::

        # From a Calibration object
        msim = sti.make_calib_sims(calib=calib, n_parsets=200)

        # From a saved parameters DataFrame
        pars_df = sc.loadobj('results/pars.df')
        msim = sti.make_calib_sims(calib_pars=pars_df, sim=make_sim(), n_parsets=50)

        # Single parameter set with multiple seeds
        msim = sti.make_calib_sims(
            calib_pars={'hiv.beta_m2f': 0.05, 'syph.beta_m2f': 0.2},
            sim=make_sim(), seeds_per_par=5,
        )

        # Different scenario with calibrated parameters
        msim = sti.make_calib_sims(
            calib_pars=pars_df, n_parsets=10, seeds_per_par=5,
            sim=make_sim(scenario='intervention'),
            check_fn=lambda s: float(np.sum(s.results.syph.new_infections[-60:])) > 0,
        )
    """
    # --- Normalize inputs to par_sets (list of flat dicts) and base_sim ---
    if calib is not None:
        par_sets = calib.get_pars(n=n_parsets)
        base_sim = sim if sim is not None else calib.sim
        if check_fn is None:
            check_fn = getattr(calib, 'check_fn', None)
    elif calib_pars is not None:
        base_sim = sim
        if isinstance(calib_pars, pd.DataFrame):
            meta = _META_KEYS
            df = calib_pars.head(n_parsets) if n_parsets else calib_pars
            par_cols = [c for c in df.columns if c not in meta]
            par_sets = [{col: row[col] for col in par_cols} for _, row in df.iterrows()]
        elif isinstance(calib_pars, dict):
            par_sets = [calib_pars]
        elif isinstance(calib_pars, list):
            par_sets = calib_pars[:n_parsets] if n_parsets else calib_pars
        else:
            raise TypeError(f'calib_pars must be a dict, list, or DataFrame, not {type(calib_pars)}')
    else:
        raise ValueError('Provide either calib or calib_pars')

    if base_sim is None:
        raise ValueError('No base sim available — pass sim= or use a Calibration that has one')

    # --- Create sims ---
    sims = []
    for i, pars in enumerate(par_sets):
        for seed in range(1, seeds_per_par + 1):
            s = sc.dcp(base_sim)
            s.pars['rand_seed'] = i * seeds_per_par + seed
            set_sim_pars(s, pars)
            s.par_idx = i
            s.seed = seed
            sims.append(s)

    # --- Run ---
    msim = ss.parallel(sims, **kwargs)

    # --- Filter ---
    if check_fn is not None:
        kept = [s for s in msim.sims if check_fn(s)]
        if seeds_per_par > 1:
            seen = set()
            filtered = []
            for s in kept:
                if s.par_idx not in seen:
                    filtered.append(s)
                    seen.add(s.par_idx)
            kept = filtered
        dropped = len(msim.sims) - len(kept)
        if dropped:
            print(f'Dropped {dropped}/{len(msim.sims)} sims via check_fn')
        msim.sims = kept

    return msim


class Calibration(ss.Calibration):
    """
    Customized STIsim calibration class.

    Inherits all the functionality of the Starsim calibration class, but adds:

    - A default build function that routes parameters using dot notation
      (e.g. ``'hiv.beta_m2f'``)
    - A default evaluation function that uses the data provided in the
      constructor
    - :meth:`get_pars` for extracting calibrated parameter sets

    If no ``build_fn`` is provided, uses :func:`default_build_fn` which looks
    up modules via ``sim.get_module()`` and sets their parameters.

    Args:
        sim (Sim): The simulation to calibrate
        calib_pars (dict): Parameters to calibrate using dot notation, e.g.
            ``{'hiv.beta_m2f': dict(low=0.01, high=0.1, guess=0.05)}``
        data (DataFrame): Calibration targets with 'time' column + result columns
        weights (dict): Optional weight multipliers per result
        extra_results (list): Additional results to track beyond data columns
        save_results (bool): Save sim results for each trial

    Examples::

        sim = make_sim()
        data = pd.read_csv('calibration_data.csv')
        calib_pars = {
            'hiv.beta_m2f': dict(low=0.01, high=0.1, guess=0.05),
            'structuredsexual.prop_f0': dict(low=0.5, high=0.9, guess=0.8),
        }
        calib = sti.Calibration(sim=sim, calib_pars=calib_pars, data=data, total_trials=100)
        calib.calibrate()

        # Extract best parameters and run multi-sim
        par_sets = calib.get_pars(n=200)
        msim = sti.make_calib_sims(calib=calib, n_parsets=200)
    """
    def __init__(self, sim, calib_pars, data=None, weights=None, extra_results=None, save_results=False, check_fn=None, **kwargs):

        # Use default build_fn if none provided
        kwargs.setdefault('build_fn', default_build_fn)

        # Flatten nested calib_pars to dot notation for Optuna
        calib_pars = flatten_calib_pars(calib_pars)

        super().__init__(sim, calib_pars, **kwargs)

        # Post-sim check function: takes a sim, returns False to reject (return inf mismatch)
        self.check_fn = check_fn

        # Custom STIsim calibration elements
        # Load data -- this is expecting a dataframe with a column for 'time' and other columns for to sim results
        self.sim_result_list = []
        self.sim_results = []
        self.weights = weights
        self.save_results = save_results
        self.results_to_save = None
        self.tmp_filename = 'tmp_calibration_%05i.obj'

        if data is not None:
            self.data = ss.validate_sim_data(data, die=True)
            self.sim_result_list = self.data.cols

            if save_results:
                self.results_to_save = list(set(self.sim_result_list + sc.tolist(extra_results)))

            # If neither components nor eval_fn have been supplied but we have a dataframe,
            # then we construct a custom eval function
            if (len(self.components) == 0) and (eval_fn not in kwargs):
                self._eval_from_data()

        return

    def _eval_from_data(self):
        """
        Make custom eval function from dataframe
        """
        result_list = self.sim_result_list if not self.save_results else self.results_to_save
        self.eval_kw = dict(data=self.data, sim_result_list=self.sim_result_list, weights=self.weights, df_res_list=result_list)
        self.eval_fn = eval_fn
        return

    def worker(self):
        """ Run a single worker, catching exceptions so one crash doesn't kill all workers """
        if self.verbose:
            op.logging.set_verbosity(op.logging.DEBUG)
        else:
            op.logging.set_verbosity(op.logging.ERROR)
        study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.study_name, sampler=self.run_args.sampler)
        try:
            output = study.optimize(self.run_trial, n_trials=self.run_args.n_trials, callbacks=None)
        except Exception as E:
            print(f'Worker failed with error: {E}')
            output = None
        return output

    def calibrate(self, calib_pars=None, **kwargs):
        """
        Perform calibration with crash-recovery support.

        If continue_db=True and the database already has completed trials, only
        the remaining trials will be run. This allows recovery from crashes by
        simply re-running the same command.
        """
        if calib_pars is not None:
            self.calib_pars = calib_pars
        self.run_args.update(kwargs)

        t0 = sc.tic()
        self.study = self.make_study()

        # Resume logic: count existing trials and only run the remainder
        if self.run_args.continue_db:
            study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.study_name)
            n_existing = len([t for t in study.trials if t.state == op.trial.TrialState.COMPLETE])
            total_trials = self.run_args.n_trials * self.run_args.n_workers
            if n_existing > 0:
                n_remaining = max(0, total_trials - n_existing)
                if n_remaining == 0:
                    print(f'Calibration already has {n_existing} completed trials, no more needed')
                else:
                    self.run_args.n_trials = int(np.ceil(n_remaining / self.run_args.n_workers))
                    print(f'Resuming calibration: {n_existing} trials complete, running ~{n_remaining} more')
                    self.run_workers()
            else:
                self.run_workers()
        else:
            self.run_workers()

        study = op.load_study(storage=self.run_args.storage, study_name=self.run_args.study_name, sampler=self.run_args.sampler)
        self.best_pars = sc.objdict(study.best_params)
        self.elapsed = sc.toc(t0, output=True)

        self.parse_study(study)

        if self.verbose:
            print('Best pars:', self.best_pars)

        self.calibrated = True
        if not self.run_args.keep_db:
            self.remove_db()

        return self

    def run_trial(self, trial):
        """ Define the objective for Optuna """
        if self.calib_pars is not None:
            pars = self._sample_from_trial(self.calib_pars, trial)
        else:
            pars = None

        if self.reseed:
            pars['rand_seed'] = trial.suggest_int('rand_seed', 0, 1_000_000) # Choose a random rand_seed

        # Prune if the prune_fn returns True
        if self.prune_fn is not None and self.prune_fn(pars):
            raise op.exceptions.TrialPruned()

        sim = self.run_sim(pars)

        # Check for invalid simulations (e.g. disease die-out)
        if sim is not None and self.check_fn is not None:
            if not self.check_fn(sim):
                return np.inf

        # Compute fit
        fit = self.eval_fn(sim, **self.eval_kw)

        if self.save_results:
            filename = self.tmp_filename % trial.number
            sc.save(filename, sim.df_res)

        return fit

    def parse_study(self, study):
        """Parse the study into a data frame -- called automatically """
        super().parse_study(study)
        if self.save_results:
            loaded = self.load_results(study)
            # Filter sim_results to only successfully loaded trials
            if loaded is not None:
                self.sim_results = [self.sim_results[i] for i in range(len(self.sim_results))]
        return

    def load_results(self, study):
        """
        Load the results from the tmp files, tracking which loaded successfully
        """
        loaded_trials = set()
        if self.save_results:
            print('Loading saved results...')
            for trial in study.trials:
                n = trial.number
                try:
                    filename = self.tmp_filename % trial.number
                    results = sc.load(filename)
                    self.sim_results.append(results)
                    loaded_trials.add(n)
                    try:
                        os.remove(filename)
                        if self.verbose: print(f'    Removed temporary file {filename}')
                    except Exception as E:
                        errormsg = f'Could not remove {filename}: {str(E)}'
                        if self.verbose: print(errormsg)
                    if self.verbose: print(f'  Loaded trial {n}')
                except Exception as E:
                    errormsg = f'Warning, could not load trial {n}: {str(E)}'
                    if self.verbose: print(errormsg)
        return loaded_trials

    def get_pars(self, n=None):
        """
        Extract top-N calibrated parameter sets as a list of flat dicts.

        Each dict maps ``'module.par'`` keys to scalar values (metadata columns
        like ``index``, ``mismatch``, and ``rand_seed`` are stripped).  The
        returned dicts can be passed directly to :func:`set_sim_pars` or
        :func:`make_calib_sims`.

        Args:
            n (int): Number of top parameter sets. ``None`` returns all.

        Returns:
            list[dict]: Parameter sets sorted by mismatch (best first).
        """
        if self.df is None:
            raise ValueError('No calibration results available — run calibrate() first')
        df = self.df.head(n) if n else self.df
        par_cols = [c for c in df.columns if c not in _META_KEYS]
        return [{col: row[col] for col in par_cols} for _, row in df.iterrows()]

    def save(self, filename, shrink=True, n_results=None, save_pars=True, pars_filename=None):
        """
        Save calibration results.

        Optionally shrinks to the top *n_results* before saving (default:
        keep top 10% of completed trials).  Also saves the parameter
        DataFrame as a separate file.

        Args:
            filename (str/Path): Path for the saved calibration object.
            shrink (bool): If True (default), shrink before saving. The full
                (unshrunk) object is saved with a ``_full`` suffix first.
            n_results (int): Number of top results to keep when shrinking.
                Default: ``len(self.df) // 10`` (minimum 1).
            save_pars (bool): If True (default), also save ``calib.df``.
            pars_filename (str/Path): Path for the parameter DataFrame. If
                None, uses ``{stem}_pars.df`` next to *filename*.

        Returns:
            The (possibly shrunk) calibration object.
        """
        filename = sc.path(filename)

        if shrink and self.sim_results:
            if n_results is None:
                n_results = max(len(self.df) // 10, 1)
            full_filename = filename.parent / (filename.stem + '_full' + filename.suffix)
            sc.save(full_filename, self)
            print(f'Saved full calibration to {full_filename}')
            calib = self.shrink(n_results=n_results)
        else:
            calib = self

        sc.save(filename, calib)
        print(f'Saved calibration to {filename}')

        if save_pars:
            if pars_filename is None:
                pars_filename = filename.parent / (filename.stem + '_pars.df')
            df = calib.df if hasattr(calib, 'df') else self.df
            sc.save(pars_filename, df)
            print(f'Saved parameters to {pars_filename}')

        return calib

    def shrink(self, n_results=100, make_df=True):
        """ Shrink the results to only the best fit """
        cal = sc.objdict()
        n_results = min(n_results, len(self.df))
        cal.sim_results = [self.sim_results[i] for i in range(n_results)]

        # Make a dataframe with the best sim and extra results
        if make_df:
            dfs = sc.autolist()
            for i in range(n_results):
                md = cal.sim_results[i]
                df = pd.DataFrame(md)
                df['res_no'] = i
                dfs += df
            cal.resdf = pd.concat(dfs)

        cal.data = self.data
        cal.df = self.df.iloc[0:n_results, ]
        return cal







