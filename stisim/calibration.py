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

__all__ = ['Calibration', 'compute_gof', 'default_build_fn']


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
    Custom evaluation function for STIsim
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


def default_build_fn(sim, calib_pars, **kwargs):
    """
    Default build function for STIsim calibration.

    Routes calibration parameters to the correct sim component based on naming
    conventions:
        - ``hiv_*``   → ``sim.diseases.hiv.pars[*]``
        - ``syph_*``  → ``sim.diseases.syphilis.pars[*]``
        - ``{disease}_*`` → ``sim.diseases.{disease}.pars[*]``
        - ``nw_*``    → ``sim.networks.structuredsexual.pars[*]``

    Parameters whose names don't match any prefix are skipped with a warning.

    Args:
        sim (Sim): An initialized simulation
        calib_pars (dict): Calibration parameters with values set by the sampler

    Returns:
        Sim: The modified simulation

    Example::

        calib_pars = dict(
            hiv_beta_m2f=dict(low=0.01, high=0.10, guess=0.035, value=0.05),
            hiv_eff_condom=dict(low=0.5, high=0.99, guess=0.95, value=0.90),
            nw_f1_conc=dict(low=0.005, high=0.3, guess=0.16, value=0.10),
        )
        sim = default_build_fn(sim, calib_pars)
    """
    if not sim.initialized:
        sim.init()

    # Build lookup of disease short names to disease objects
    disease_map = {}
    for name, disease in sim.diseases.items():
        disease_map[name] = disease
        # Also map common short names
        if name == 'syphilis':
            disease_map['syph'] = disease

    # Get the structured sexual network if available
    nw = None
    if 'structuredsexual' in sim.networks:
        nw = sim.networks.structuredsexual

    for k, pars in calib_pars.items():
        # Skip metadata keys
        if k in ['rand_seed', 'index', 'mismatch']:
            continue

        # Extract value
        if isinstance(pars, dict):
            v = pars.get('value', pars)
        elif sc.isnumber(pars):
            v = pars
        else:
            continue

        # Route by prefix
        matched = False

        # Check disease prefixes
        for prefix, disease in disease_map.items():
            if k.startswith(prefix + '_'):
                par_name = k[len(prefix) + 1:]
                disease.pars[par_name] = v
                matched = True
                break

        # Check network prefix
        if not matched and k.startswith('nw_') and nw is not None:
            par_name = k[3:]
            if hasattr(nw.pars[par_name], 'set'):
                nw.pars[par_name].set(v)
            else:
                nw.pars[par_name] = v
            matched = True

        if not matched:
            print(f'Warning: calibration parameter "{k}" not matched to any disease or network')

    return sim


class Calibration(ss.Calibration):
    """
    Customized STIsim calibration class.

    Inherits all the functionality of the Starsim calibration class, but adds:
    - A default build function that routes parameters by prefix (e.g. ``hiv_beta_m2f``)
    - A default evaluation function that uses the data provided in the constructor

    If no ``build_fn`` is provided, uses :func:`default_build_fn` which routes
    parameters to diseases and networks by prefix convention.

    Args:
        sim (Sim): The simulation to calibrate
        calib_pars (dict): Parameters to calibrate, e.g. ``dict(hiv_beta_m2f=dict(low=0.01, high=0.1))``
        data (DataFrame): Calibration targets with 'time' column + result columns
        weights (dict): Optional weight multipliers per result
        extra_results (list): Additional results to track beyond data columns
        save_results (bool): Save sim results for each trial

    Examples::

        sim = make_sim()
        data = pd.read_csv('calibration_data.csv')
        calib_pars = dict(
            hiv_beta_m2f=dict(low=0.01, high=0.1, guess=0.05),
            nw_prop_f0=dict(low=0.5, high=0.9, guess=0.8),
        )
        calib = stisim.Calibration(sim=sim, calib_pars=calib_pars, data=data, total_trials=100)
        calib.calibrate()
    """
    def __init__(self, sim, calib_pars, data=None, weights=None, extra_results=None, save_results=False, check_fn=None, **kwargs):

        # Use default build_fn if none provided
        kwargs.setdefault('build_fn', default_build_fn)

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







