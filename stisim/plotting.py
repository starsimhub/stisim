"""
Default plotting functions for STIsim.
"""
import numpy as np
import sciris as sc
import pylab as pl

__all__ = ['plot_hiv', 'plot_hiv_msim']

# Default result keys for HIV plotting
HIV_PLOT_KEYS = [
    'hiv_new_infections',
    'hiv_new_deaths',
    'hiv_n_infected',
    'hiv_prevalence',
    'hiv_n_on_art',
    'n_alive',
]

HIV_PLOT_TITLES = [
    'New HIV infections',
    'HIV-related deaths',
    'PLHIV',
    'HIV prevalence (%)',
    'Number on ART',
    'Population size',
]

# Mapping from result key (underscore) to data column name.
# After validate_sim_data, the time column becomes the index,
# and the data columns may use dot or underscore notation.
_DATA_COL_MAP = {
    'hiv_new_infections': ['hiv.new_infections', 'hiv_new_infections'],
    'hiv_new_deaths':     ['hiv.new_deaths', 'hiv_new_deaths'],
    'hiv_n_infected':     ['hiv.n_infected', 'hiv_n_infected'],
    'hiv_prevalence':     ['hiv.prevalence', 'hiv_prevalence'],
    'hiv_n_on_art':       ['hiv.n_on_art', 'hiv_n_on_art'],
    'n_alive':            ['n_alive'],
}


def _find_data_col(data, key):
    """Find a matching column in data for a given result key."""
    candidates = _DATA_COL_MAP.get(key, [key])
    for col in candidates:
        if col in data.columns:
            return col
    return None


def plot_hiv(sim, keys=None, start_year=None, end_year=None, data=None, fig_kw=None,
             do_save=False, filename=None, folder=None, dpi=150):
    """
    Plot key HIV simulation results in a 2x3 panel grid.

    Overlays comparison data (e.g. UNAIDS estimates) if available via sim.data.

    Args:
        sim (Sim): A completed stisim/hivsim simulation
        keys (list): Result keys to plot (default: standard HIV panel)
        start_year (float): First year to show (default: sim start)
        end_year (float): Last year to show (default: sim end)
        data (DataFrame): Comparison data to overlay; if None, uses sim.data
        fig_kw (dict): Passed to pl.subplots()
        do_save (bool): Whether to save the figure to disk
        filename (str): Filename for saving (default: 'hiv_plot.png')
        folder (str): Folder for saving (default: current directory)
        dpi (int): Resolution for saved figure

    Returns:
        matplotlib.figure.Figure

    Examples:
        >>> import stisim as sti
        >>> sim = sti.Sim(diseases='hiv', n_agents=5000, dur=30)
        >>> sim.run()
        >>> sti.plot_hiv(sim)
        >>> sti.plot_hiv(sim, do_save=True, filename='my_plot.png', folder='figures')
    """
    if not sim.results_ready:
        raise RuntimeError('Please run the sim before plotting')

    if keys is None:
        keys = HIV_PLOT_KEYS
    titles = {k: t for k, t in zip(HIV_PLOT_KEYS, HIV_PLOT_TITLES)}

    if data is None:
        data = getattr(sim, 'data', None)

    # Get the results as a dataframe
    df = sim.to_df(resample='year', use_years=True, sep='_')
    if start_year is None:
        start_year = df['timevec'].min()
    if end_year is None:
        end_year = df['timevec'].max()
    dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
    dfplot = dfplot.set_index('timevec')

    # Filter data to time range (after validate_sim_data, index is the time column)
    if data is not None:
        data = data.loc[(data.index >= start_year) & (data.index <= end_year)]

    # Create figure
    n_keys = len(keys)
    nrows, ncols = sc.getrowscols(n_keys)
    fig_defaults = dict(figsize=(6*ncols, 4*nrows))
    fig_kw = sc.mergedicts(fig_defaults, fig_kw)
    fig, axes = pl.subplots(nrows, ncols, **fig_kw)
    axes = np.array(axes).ravel()

    for i, key in enumerate(keys):
        ax = axes[i]

        # Plot data overlay
        if data is not None:
            data_col = _find_data_col(data, key)
            if data_col is not None:
                vals = data[data_col].dropna()
                if len(vals) > 0:
                    yvals = vals * 100 if key == 'hiv_prevalence' else vals
                    ax.scatter(vals.index, yvals, color='k', zorder=10, s=15, label='Data')

        # Plot sim results
        if key in dfplot.columns:
            x = dfplot.index
            y = dfplot[key]
            if key == 'hiv_prevalence':
                y = y * 100
            ax.plot(x, y, label='Model')

        ax.set_title(titles.get(key, key))
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)

    # Clean up unused axes
    for j in range(n_keys, len(axes)):
        axes[j].set_visible(False)

    sc.figlayout(fig=fig)

    if do_save:
        if filename is None:
            filename = 'hiv_plot.png'
        if folder is not None:
            filename = sc.path(folder) / filename
        sc.savefig(str(filename), fig=fig, dpi=dpi)

    return fig


def plot_hiv_msim(msim, keys=None, start_year=None, end_year=None, data=None,
                  quantiles=None, fig_kw=None, do_save=False, filename=None,
                  folder=None, dpi=150):
    """
    Plot MultiSim HIV results with median line and shaded IQR.

    Args:
        msim (MultiSim): A completed MultiSim
        keys (list): Result keys to plot (default: standard HIV panel)
        start_year (float): First year to show (default: sim start)
        end_year (float): Last year to show (default: sim end)
        data (DataFrame): Comparison data to overlay; if None, uses base_sim.data
        quantiles (dict): Quantile bounds, e.g. {'low': 0.25, 'high': 0.75}
        fig_kw (dict): Passed to pl.subplots()
        do_save (bool): Whether to save the figure to disk
        filename (str): Filename for saving (default: 'hiv_msim_plot.png')
        folder (str): Folder for saving (default: current directory)
        dpi (int): Resolution for saved figure

    Returns:
        matplotlib.figure.Figure
    """
    import pandas as pd

    if not msim.sims or not msim.sims[0].results_ready:
        raise RuntimeError('Please run the MultiSim before plotting')

    if keys is None:
        keys = HIV_PLOT_KEYS
    titles = {k: t for k, t in zip(HIV_PLOT_KEYS, HIV_PLOT_TITLES)}

    if quantiles is None:
        quantiles = {'low': 0.25, 'high': 0.75}

    if data is None:
        data = getattr(msim.base_sim, 'data', None)

    # Collect dataframes from all sims
    dfs = []
    for sim in msim.sims:
        df = sim.to_df(resample='year', use_years=True, sep='_')
        dfs.append(df)

    # Compute median and quantiles
    combined = pd.concat(dfs)
    grouped = combined.groupby('timevec')
    df_median = grouped.median(numeric_only=True)
    df_lo = grouped.quantile(quantiles['low'], numeric_only=True)
    df_hi = grouped.quantile(quantiles['high'], numeric_only=True)

    if start_year is None:
        start_year = df_median.index.min()
    if end_year is None:
        end_year = df_median.index.max()
    df_median = df_median.loc[(df_median.index >= start_year) & (df_median.index <= end_year)]
    df_lo = df_lo.loc[df_median.index]
    df_hi = df_hi.loc[df_median.index]

    # Filter data
    if data is not None:
        data = data.loc[(data.index >= start_year) & (data.index <= end_year)]

    # Create figure
    n_keys = len(keys)
    nrows, ncols = sc.getrowscols(n_keys)
    fig_defaults = dict(figsize=(6*ncols, 4*nrows))
    fig_kw = sc.mergedicts(fig_defaults, fig_kw)
    fig, axes = pl.subplots(nrows, ncols, **fig_kw)
    axes = np.array(axes).ravel()

    for i, key in enumerate(keys):
        ax = axes[i]
        x = df_median.index
        scale = 100 if key == 'hiv_prevalence' else 1

        # Data overlay
        if data is not None:
            data_col = _find_data_col(data, key)
            if data_col is not None:
                vals = data[data_col].dropna()
                if len(vals) > 0:
                    ax.scatter(vals.index, vals * scale, color='k', zorder=10, s=15, label='Data')

        # Median + IQR
        if key in df_median.columns:
            y_med = df_median[key] * scale
            y_lo = df_lo[key] * scale
            y_hi = df_hi[key] * scale
            ax.fill_between(x, y_lo, y_hi, alpha=0.3, label=f'{int(quantiles["low"]*100)}-{int(quantiles["high"]*100)}%')
            ax.plot(x, y_med, label='Median')

        ax.set_title(titles.get(key, key))
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)

    for j in range(n_keys, len(axes)):
        axes[j].set_visible(False)

    sc.figlayout(fig=fig)

    if do_save:
        if filename is None:
            filename = 'hiv_msim_plot.png'
        if folder is not None:
            filename = sc.path(folder) / filename
        sc.savefig(str(filename), fig=fig, dpi=dpi)

    return fig
