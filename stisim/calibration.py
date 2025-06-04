"""
Define the calibration class
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import os


__all__ = ['Calibration', 'compute_gof']


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
                rescomps = sres.split('_')
                modname = rescomps[0]
                if isinstance(df_sim.results[modname], ss.Results):
                    resname = '_'.join(rescomps[1:])
                    dfs += df_sim.results[modname][resname].to_df(resample='year', use_years=True, col_names=sres)

    df_res = pd.concat(dfs, axis=1)
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


class Calibration(ss.Calibration):
    """
    Customized STIsim calibration class
    Inherits all the functionality of the Starsim calibration class, but adds a
    default evaluation function that uses the data provided in the constructor
    """
    def __init__(self, sim, calib_pars, data=None, weights=None, extra_results=None, save_results=False, **kwargs):

        super().__init__(sim, calib_pars, **kwargs)

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

        # Compute fit
        fit = self.eval_fn(sim, **self.eval_kw)

        if self.save_results:
            filename = self.tmp_filename % trial.number
            sc.save(filename, sim.df_res)

        return fit

    def parse_study(self, study):
        """Parse the study into a data frame -- called automatically """
        super().parse_study(study)
        self.load_results(study)
        self.sim_results = [self.sim_results[i] for i in self.df.index]
        return

    def load_results(self, study):
        """
        Load the results from the tmp files
        """
        if self.save_results:
            print('Loading saved results...')
            for trial in study.trials:
                n = trial.number
                try:
                    filename = self.tmp_filename % trial.number
                    results = sc.load(filename)
                    self.sim_results.append(results)
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
        return

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
