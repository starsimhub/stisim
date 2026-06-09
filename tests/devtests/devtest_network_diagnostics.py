"""
Network dynamics validation: verify that network parameters (concurrency, pair
formation, relationship duration, debut age, MSM) affect behaviour as expected.
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sciris as sc
import shutil
import starsim as ss
import stisim as sti
import sys
import time

from stisim.analyzers import PartnershipFormationAnalyzer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from stisim.networks.matchers import MATCHERS

# All output files written by these diagnostics go in this single directory.
OUTPUT_DIR = sc.thisdir(__file__, 'devtest_network_diagnostics')


def reset_output_dir():
    """Delete OUTPUT_DIR if it exists and recreate it empty."""
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    return


def _output_path(filename):
    """Return the full path for an output file inside OUTPUT_DIR."""
    return os.path.join(OUTPUT_DIR, filename)


# Reset the output directory once when the module is imported, so a clean
# directory is guaranteed whether tests run via pytest or via __main__.
reset_output_dir()

@sc.timer()
def test_age_differences(n_agents=25000, target_age_gap=8, sd=3, matching_algo=MATCHERS['closest_age_tapered_seeking']):
    """
    Sexual network with default relationship type distribution: run for 1 month and
    analyze partner age differences. Produces a scatterplot of female vs male partner
    ages with a best-fit regression line, and prints the mean and standard deviation
    of the age gap.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt


    age_diff_pars_older = {'teens': [(target_age_gap, sd), (target_age_gap, sd), (target_age_gap, sd)],
                           'young': [(target_age_gap, sd), (target_age_gap, sd), (target_age_gap, sd)],
                           'adult': [(target_age_gap, sd), (target_age_gap, sd), (target_age_gap, sd)]}
    network = sti.MFNetwork(age_diff_pars=age_diff_pars_older, match_method=matching_algo)
    sim = sti.Sim(n_agents=n_agents, networks=[network], dur=1, rand_seed=1)
    sim.run()

    nw = sim.networks.mfnetwork
    male_ages = np.array(nw.edges.age_p1, dtype=float)
    female_ages = np.array(nw.edges.age_p2, dtype=float)

    # Restrict both axes to 15–75 years
    min_age = 0
    max_age = 100
    mask = (female_ages >= min_age) & (female_ages <= max_age) & (male_ages >= min_age) & (male_ages <= max_age)
    female_ages = female_ages[mask]
    male_ages = male_ages[mask]

    assert len(female_ages) > 0, "No partnerships found after 1 month — check network setup"

    age_diffs = male_ages - female_ages
    mean_diff = np.mean(age_diffs)
    std_diff = np.std(age_diffs)
    print(f"Average partner age difference (male − female): {mean_diff:.2f} years")
    print(f"Standard deviation of partner age difference:   {std_diff:.2f} years")

    # Linear regression via numpy
    coeffs = np.polyfit(female_ages, male_ages, 1)
    slope, intercept = coeffs

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(female_ages, male_ages, alpha=0.4, s=8, label='Partnerships')
    x_fit = np.array([min_age, max_age])
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'r-', linewidth=2,
            label=f'Linear fit  slope={slope:.2f}')
    ax.set_xlabel('Female partner age (years)')
    ax.set_ylabel('Male partner age (years)')
    ax.set_title('Relationship partner ages (at formation time)')
    ax.set_xlim(min_age, max_age)
    ax.set_ylim(min_age, max_age)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outfile = _output_path('relationship_age_differences.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Scatterplot saved to {outfile}")

    return sim



def _run_age_diff_sim(target_age_gap, seed, algorithm, n_agents=25000, dur=1):
    age_diff_pars = {
        'teens': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
        'young': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
        'adult': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
    }
    network = sti.MFNetwork(age_diff_pars=age_diff_pars, match_method=algorithm)
    # with warnings.catch_warnings():
    #     warnings.simplefilter('ignore')
    sim = sti.Sim(n_agents=n_agents, networks=[network], dur=dur,
                  rand_seed=seed) # , verbose=0)
    t_start = time.time()
    sim.run()
    t_end = time.time()
    run_time = t_end - t_start

    nw = sim.networks.mfnetwork
    male_ages = np.array(nw.edges.age_p1, dtype=float)
    female_ages = np.array(nw.edges.age_p2, dtype=float)

    if len(male_ages) == 0:
        return float('nan'), float('nan'), 0
    age_diffs = male_ages - female_ages
    return float(np.mean(age_diffs)), float(np.std(age_diffs)), int(len(age_diffs)), run_time


def _run_age_diff_sweep(target_age_gaps, seeds, algorithm, n_agents=10000, dur=1, label=''):
    """Run an (upper_age x seed) grid and return a results dict.

    Prints one line per (upper_age, seed), returns:
        {upper_age: {'means': [..], 'stds': [..], 'n_pairs': [..]}}
    """
    print(f"  [{label}] n_agents={n_agents} dur={dur}")
    results = {}
    for target_age_gap in target_age_gaps:
        means, stds, ns, run_times = [], [], [], []
        for seed in seeds:
            mean, std, n, run_time = _run_age_diff_sim(target_age_gap, seed, algorithm, n_agents=n_agents, dur=dur)
            means.append(mean)
            stds.append(std)
            ns.append(n)
            run_times.append(run_time)
            print(f"  [{label}] target_gap={target_age_gap:2d} seed={seed}: "
                  f"n={n:5d} mean={mean:7.3f} std={std:6.3f}")
        results[target_age_gap] = {'means': means, 'stds': stds, 'n_pairs': ns, 'time_mean': np.mean(run_times), 'time_stds': np.std(run_times)}
    return results


def _print_sweep_summary(results, label=''):
    """Print a per-upper_age summary: mean over seeds of mean and of std."""
    print(f"\n{label} summary (mean over seeds):")
    print(f"  target_age | mean(age_diff) | mean(std) | mean(n_pairs)")
    print(f"  ----------+----------------+-----------+--------------")
    for target_gap, vals in results.items():
        mm = float(np.nanmean(vals['means']))
        ms = float(np.nanmean(vals['stds']))
        mn = int(np.mean(vals['n_pairs']))
        print(f"  {target_gap:9d} | {mm:14.3f} | {ms:9.3f} | {mn:13d}")


@sc.timer()
def test_algorithm_comparison(n_agents=25000):
    """this test is useful for comparing new proposed algorithms to the current one"""
    target_age_gaps = [5, 8, 11]
    seeds = [1, 2, 3, 4, 5]
    dur = 25

    algorithms = {'closest_age_tapered_seeking': MATCHERS['closest_age_tapered_seeking']}

    all_results = {}
    for matching_algo_name, matching_algo in algorithms.items():
        print(f"\n=== Sweep: {matching_algo} ===")
        all_results[matching_algo_name] = _run_age_diff_sweep(target_age_gaps, seeds, algorithm=matching_algo,
                                                     n_agents=n_agents, dur=dur, label=matching_algo_name)
        _print_sweep_summary(all_results[matching_algo_name], label=matching_algo_name)

    # Side-by-side: mean age diff (mean over seeds), one column per algorithm
    print("\n=== Side-by-side: mean age diff per target_gap ===")
    header = "  target_gap | " + " | ".join(f"{name:>30s}" for name in algorithms.keys())
    sep    = "  " + "-"*9 + "-+-" + "-+-".join(["-"*30]*len(algorithms.keys()))
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(np.nanmean(all_results[name][target_age_gap]['means'])):30.3f}"
                 for name in algorithms.keys()]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean realized std per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(np.nanmean(all_results[name][target_age_gap]['stds'])):30.3f}"
                 for name in algorithms.keys()]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean partnership count per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{int(np.mean(all_results[name][target_age_gap]['n_pairs'])):30d}"
                 for name in algorithms.keys()]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean run time (s) per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(all_results[name][target_age_gap]['time_mean']):30.3f}"
                 for name in algorithms.keys()]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: run time std (s) per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(all_results[name][target_age_gap]['time_stds']):30.3f}"
                 for name in algorithms.keys()]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    return {'stats': all_results}


@sc.timer()
def test_n_partners_distribution(n_agents=25000, n_runs=5, dur=25, window_months=12, target_age_gap=8,
                                 matching_algo=MATCHERS['closest_age_tapered_seeking']):
    """
    Sweep over rand_seed and characterize the distribution of unique female
    partners per male over a trailing 12-month window in an MFNetwork.

    Produces a histogram whose bar heights are the mean normalized fraction
    of males per partner-count bin (across the seed sweep), with 5%/95%
    normal-distribution whiskers, and a table of equivalent data including
    the underlying mean male counts per bin.

    Args:
        n_agents (int): Number of agents per simulation.
        n_runs (int): Number of seeds to sweep over.
        dur (float): Simulation duration in years.
        window_months (int): Trailing partner-count window in months.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    sdA = 3
    age_diff_pars = {'teens': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'young': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'adult': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)]}

    seeds = list(range(n_runs))
    n_explicit_bins = 10  # bins 0..9 explicit, plus a catchall "10+"
    bin_edges = np.concatenate([np.arange(n_explicit_bins + 1), [np.inf]])
    n_bins = n_explicit_bins + 1
    bin_labels = [str(i) for i in range(n_explicit_bins)] + [f'{n_explicit_bins}+']

    per_run_counts = np.zeros((n_runs, n_bins), dtype=float)
    per_run_n_males = np.zeros(n_runs, dtype=int)

    for i, seed in enumerate(seeds):
        network = sti.MFNetwork(age_diff_pars=age_diff_pars, match_method=matching_algo)
        analyzer = PartnershipFormationAnalyzer(networks=['mfnetwork'])
        sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                      dur=dur, rand_seed=seed)
        sim.run()
        # sti.Sim copies the analyzer during init, so retrieve the live instance
        sim_analyzer = sim.analyzers.partnershipformationanalyzer
        counts = sim_analyzer.get_unique_partners_active_per_agent(female=False, window_months=window_months)['mfnetwork']
        binned, _ = np.histogram(counts, bins=bin_edges)
        per_run_counts[i] = binned
        per_run_n_males[i] = counts.size
        print(f"  seed={seed}: n_males={counts.size} "
              f"mean_partners={counts.mean():.3f} max_partners={int(counts.max()) if counts.size else 0}")

    # Per-run normalized fractions, then aggregate across runs
    safe_n_males = np.where(per_run_n_males > 0, per_run_n_males, 1)
    per_run_fractions = per_run_counts / safe_n_males[:, None]
    mean_fractions = per_run_fractions.mean(axis=0)
    std_fractions = per_run_fractions.std(axis=0)
    ci_low = mean_fractions - 1.645 * std_fractions
    ci_high = mean_fractions + 1.645 * std_fractions

    # Underlying raw-count stats for the table
    mean_counts = per_run_counts.mean(axis=0)
    std_counts = per_run_counts.std(axis=0)

    # Histogram plot (dependent axis: normalized fraction of males)
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(n_bins)
    ax.bar(x, mean_fractions, color='steelblue', edgecolor='black', alpha=0.8)
    yerr = np.vstack([mean_fractions - ci_low, ci_high - mean_fractions])
    yerr = np.clip(yerr, 0, None)
    ax.errorbar(x, mean_fractions, yerr=yerr, fmt='none', ecolor='black',
                capsize=4, capthick=1.5, lw=1.5)
    for xi, mf in zip(x, mean_fractions):
        if mf > 0:
            ax.text(xi, mf / 2, f'{mf:.3f}', ha='center', va='bottom',
                    fontsize=8, color='black')
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel(f'Number of unique female partners (trailing {window_months} months)')
    ax.set_ylabel('Mean fraction of males (across seed sweep)')
    ax.set_title(f'Partner-count distribution per male '
                 f'({n_runs} runs, {n_agents} agents, {dur} yr)')
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    outfile = _output_path('n_partners_distribution.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Histogram saved to {outfile}")

    print(f"\n  partners | mean(frac) | std(frac) |  5% CI   |  95% CI  | mean(n_males) | std(n_males)")
    print(f"  ---------+------------+-----------+----------+----------+---------------+--------------")
    for label, mf, sf, lo, hi, mc, sd in zip(bin_labels, mean_fractions, std_fractions, ci_low, ci_high, mean_counts, std_counts):
        print(f"  {label:>8s} | {mf:10.4f} | {sf:9.4f} | {lo:8.4f} | {hi:8.4f} | {mc:13.2f} | {sd:12.2f}")

    return {'mean_fractions': mean_fractions, 'std_fractions': std_fractions,
            'ci_low': ci_low, 'ci_high': ci_high,
            'mean_counts': mean_counts, 'std_counts': std_counts,
            'bin_labels': bin_labels,
            'per_run_counts': per_run_counts, 'per_run_fractions': per_run_fractions}


@sc.timer()
def test_n_partners_pyramid(n_agents=25000, n_runs=5, dur=25, window_months=12, target_age_gap=8,
                            matching_algo=MATCHERS['closest_age_tapered_seeking']):
    """
    Population-pyramid of partnerships formed (non-unique) by age and sex.

    Sweeps over ``n_runs`` seeds and, for each, uses ``PartnershipFormationAnalyzer``
    to count the number of partnerships *formed* (not unique-partner counts) in the
    trailing ``window_months``, binned by each subject's age at formation. The
    per-bin counts are averaged across seeds and drawn as a population pyramid:
    a central vertical line at x=0 with male counts extending left and female
    counts extending right, age bins stacked bottom (youngest) to top (oldest).
    Each bar is labelled with its raw (mean) count.

    Args:
        n_agents (int): Number of agents per simulation.
        n_runs (int): Number of seeds to sweep over.
        dur (float): Simulation duration in years.
        window_months (int): Trailing formation window in months.
        target_age_gap (float): Mean male-female age gap for partner matching.
        matching_algo: Partner-matching algorithm (from ``MATCHERS``).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    sdA = 3
    age_diff_pars = {'teens': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'young': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'adult': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)]}

    # 5-year age bins, low -> high (bottom -> top of the pyramid).
    age_bins = [f'{lo}-{lo + 5}' for lo in range(0, 80, 5)]  # '0-5' .. '75-80'
    n_bins = len(age_bins)

    seeds = list(range(n_runs))
    male_per_run = np.zeros((n_runs, n_bins), dtype=float)
    female_per_run = np.zeros((n_runs, n_bins), dtype=float)

    for i, seed in enumerate(seeds):
        network = sti.MFNetwork(age_diff_pars=age_diff_pars, match_method=matching_algo)
        analyzer = PartnershipFormationAnalyzer(networks=['mfnetwork'])
        sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                      dur=dur, rand_seed=seed)
        sim.run()
        # sti.Sim copies the analyzer during init, so retrieve the live instance
        sim_analyzer = sim.analyzers.partnershipformationanalyzer
        # {nw: {age_bin: array([window_sum])}}; collapse each bin's length-1 array.
        male_by_bin = sim_analyzer.get_n_partnerships_formed(
            female=False, age_bins=age_bins, window_months=window_months)['mfnetwork']
        female_by_bin = sim_analyzer.get_n_partnerships_formed(
            female=True, age_bins=age_bins, window_months=window_months)['mfnetwork']
        male_per_run[i] = [male_by_bin[b][0] for b in age_bins]
        female_per_run[i] = [female_by_bin[b][0] for b in age_bins]
        print(f"  seed={seed}: male_formed={int(male_per_run[i].sum())} "
              f"female_formed={int(female_per_run[i].sum())}")

    male_mean = male_per_run.mean(axis=0)
    female_mean = female_per_run.mean(axis=0)
    # Population std across the n_runs seeds (per age bin, per sex).
    male_std = male_per_run.std(axis=0)
    female_std = female_per_run.std(axis=0)
    # Two-tailed 5%/95% whiskers under a normal approximation: mean +/- z*std,
    # z=1.645 puts the lower cap at the 5th percentile and the upper cap at the
    # 95th percentile (a 90% central interval). Matches test_n_partners_distribution.
    Z_90 = 1.645
    male_err = Z_90 * male_std
    female_err = Z_90 * female_std

    # Population pyramid: males extend left (negative), females extend right.
    y = np.arange(n_bins)
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.barh(y, -male_mean, color='blue', alpha=0.7, edgecolor='black', label='Male')
    ax.barh(y, female_mean, color='red', alpha=0.7, edgecolor='black', label='Female')
    # 5%/95% normal-approx whiskers at each bar tip (computed across seeds).
    ax.errorbar(-male_mean, y, xerr=male_err, fmt='none', ecolor='black',
                capsize=3, capthick=1.0, lw=1.0)
    ax.errorbar(female_mean, y, xerr=female_err, fmt='none', ecolor='black',
                capsize=3, capthick=1.0, lw=1.0)
    ax.axvline(0, color='black', linewidth=1.0)

    # Raw (mean) count displayed within each box.
    for yi, (m, f) in enumerate(zip(male_mean, female_mean)):
        if m > 0:
            ax.text(-m / 2, yi, f'{m:.1f}', ha='center', va='center', fontsize=8, color='black')
        if f > 0:
            ax.text(f / 2, yi, f'{f:.1f}', ha='center', va='center', fontsize=8, color='black')

    # Symmetric x-limits: abs(left edge) == right edge, large enough to contain
    # every drawn element (bar tip + its 95% whisker) on both sides.
    max_extent = max((male_mean + male_err).max(), (female_mean + female_err).max())
    lim = max_extent * 1.05 if max_extent > 0 else 1.0  # 5% padding
    ax.set_xlim(-lim, lim)

    ax.set_yticks(y)
    ax.set_yticklabels(age_bins)
    ax.set_ylabel('Age bin at formation')
    ax.set_xlabel('Mean partnerships formed (male ← | → female)')
    # x tick labels as absolute magnitudes (both directions read positive).
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f'{abs(v):g}'))
    matcher_name = getattr(matching_algo, '__name__', str(matching_algo))
    ax.set_title(f'Partnerships formed by age & sex (population pyramid)\n'
                 f'n_agents={n_agents}, n_runs={n_runs}, dur={dur}, '
                 f'window_months={window_months}, target_age_gap={target_age_gap}, '
                 f'matcher={matcher_name}')
    ax.legend(loc='upper right')
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    outfile = _output_path('n_partners_pyramid.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Pyramid saved to {outfile}")

    return {'age_bins': age_bins, 'male_mean': male_mean, 'female_mean': female_mean,
            'male_std': male_std, 'female_std': female_std,
            'male_per_run': male_per_run, 'female_per_run': female_per_run}


class AgeGapAnalyzer(ss.Analyzer):
    """Track male-female age gaps of unique relationships over a trailing window.

    Walks the active MFNetwork edges each step and records, for each unique
    (male, female) pair seen in the trailing ``window_months`` months, the
    male-minus-female age gap captured at the first appearance of that edge.
    After the sim finishes, ``get_age_gaps`` returns the per-relationship gaps.

    Args:
        network (str): Network name to inspect (default ``'mfnetwork'``).
        window_months (int): Trailing window in months; assumes monthly dt
            (one ti per month) so this is also the number of trailing steps.
    """

    def __init__(self, network='mfnetwork', window_months=12, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.network = network
        self.window_months = window_months
        # list of (ti, p1_uids, p2_uids, age_p1, age_p2)
        self._records = []
        self._final_ti = None
        return

    def step(self):
        nw = self.sim.networks[self.network]
        p1 = np.asarray(nw.edges.p1, dtype=np.int64)
        p2 = np.asarray(nw.edges.p2, dtype=np.int64)
        age_p1 = np.asarray(nw.edges.age_p1, dtype=float)
        age_p2 = np.asarray(nw.edges.age_p2, dtype=float)
        if len(p1):
            self._records.append((self.ti, p1, p2, age_p1, age_p2))
        self._final_ti = self.ti
        return

    def get_age_gaps(self):
        """Return per-relationship male-female age gaps over the trailing window."""
        threshold = self._final_ti - self.window_months + 1
        seen = {}
        for ti, p1s, p2s, a1s, a2s in self._records:
            if ti < threshold:
                continue
            for m, f, am, af in zip(p1s, p2s, a1s, a2s):
                pair = (int(m), int(f))
                if pair not in seen:
                    seen[pair] = float(am) - float(af)
        return np.array(list(seen.values()), dtype=float)


@sc.timer()
def test_age_gap_distribution(n_agents=25000, n_runs=5, dur=25, window_months=12,
                              target_age_gap=8, gap_bin_width=1, matching_algo=MATCHERS['closest_age_tapered_seeking']):
    """
    Sweep over rand_seed and characterize the distribution of male-minus-female
    age gaps across unique relationships over a trailing 12-month window in an
    MFNetwork.

    Produces a histogram whose bar heights are the mean normalized fraction of
    relationships per age-gap bin (across the seed sweep), with 5%/95%
    normal-distribution whiskers, and a table of equivalent data including the
    underlying mean relationship counts per bin. Bins are integer-centered and
    trimmed to ``[-gap_range, +gap_range]``; relationships with gaps outside
    that range contribute to the denominator but are not shown.

    Args:
        n_agents (int): Number of agents per simulation.
        n_runs (int): Number of seeds to sweep over.
        dur (float): Simulation duration in years.
        window_months (int): Trailing relationship window in months.
        gap_range (float): Half-width of the displayed gap range, in years.
        gap_bin_width (float): Histogram bin width in years.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    sdA = 3
    age_diff_pars = {'teens': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'young': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)],
                     'adult': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdA)]}
    max_allowed_age_delta = target_age_gap + 3 * sdA


    seeds = list(range(1, n_runs + 1))
    # Bins of width ``gap_bin_width`` centered on integers from
    # -max_allowed_age_delta to +max_allowed_age_delta inclusive. With
    # gap_bin_width=1 and max_allowed_age_delta=11 that yields 23 bins
    # centered on -11, -10, ..., +11 with edges at -11.5, -10.5, ..., +11.5.
    n_bins = int(round(2 * max_allowed_age_delta / gap_bin_width)) + 1
    half = gap_bin_width / 2
    bin_edges = np.linspace(-max_allowed_age_delta - half,
                             max_allowed_age_delta + half,
                             n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    if gap_bin_width == 1:
        bin_labels = [f'{int(round(c))}' for c in bin_centers]
    else:
        bin_labels = [f'{c:.1f}' for c in bin_centers]

    per_run_counts = np.zeros((n_runs, n_bins), dtype=float)
    per_run_n_rel = np.zeros(n_runs, dtype=int)
    per_run_gap_mean = np.full(n_runs, np.nan)
    per_run_gap_std = np.full(n_runs, np.nan)
    all_gaps = []

    for i, seed in enumerate(seeds):
        network = sti.MFNetwork(age_diff_pars=age_diff_pars, match_method=matching_algo)

        analyzer = AgeGapAnalyzer(network='mfnetwork', window_months=window_months)
        sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                      dur=dur, rand_seed=seed)
        sim.run()
        # sti.Sim copies the analyzer during init, so retrieve the live instance
        sim_analyzer = sim.analyzers.agegapanalyzer
        gaps = sim_analyzer.get_age_gaps()
        binned, _ = np.histogram(gaps, bins=bin_edges)
        per_run_counts[i] = binned
        per_run_n_rel[i] = gaps.size
        mean_gap = float(gaps.mean()) if gaps.size else float('nan')
        std_gap = float(gaps.std()) if gaps.size else float('nan')
        per_run_gap_mean[i] = mean_gap
        per_run_gap_std[i] = std_gap
        if gaps.size:
            all_gaps.append(gaps)
        print(f"  seed={seed}: n_rel={gaps.size} mean_gap={mean_gap:.3f} std_gap={std_gap:.3f}")

    # Per-run normalized fractions, then aggregate across runs
    safe_n_rel = np.where(per_run_n_rel > 0, per_run_n_rel, 1)
    per_run_fractions = per_run_counts / safe_n_rel[:, None]
    mean_fractions = per_run_fractions.mean(axis=0)
    std_fractions = per_run_fractions.std(axis=0)
    ci_low = mean_fractions - 1.645 * std_fractions
    ci_high = mean_fractions + 1.645 * std_fractions

    mean_counts = per_run_counts.mean(axis=0)
    std_counts = per_run_counts.std(axis=0)

    # Histogram plot (dependent axis: normalized fraction of relationships)
    fig, ax = plt.subplots(figsize=(13, 6))
    x = np.arange(n_bins)
    ax.bar(x, mean_fractions, color='steelblue', edgecolor='black', alpha=0.8)
    yerr = np.vstack([mean_fractions - ci_low, ci_high - mean_fractions])
    yerr = np.clip(yerr, 0, None)
    ax.errorbar(x, mean_fractions, yerr=yerr, fmt='none', ecolor='black',
                capsize=4, capthick=1.5, lw=1.5)
    for xi, mf in zip(x, mean_fractions):
        if mf > 0:
            ax.text(xi, mf / 2, f'{mf:.3f}', ha='center', va='bottom',
                    fontsize=7, color='black', rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=0 if n_bins <= 16 else 45)
    ax.set_xlabel(f'Male - female age gap, years (trailing {window_months} months)')
    ax.set_ylabel('Mean fraction of relationships (across seed sweep)')
    ax.set_title(f'Age-gap distribution across relationships '
                 f'({n_runs} runs, {n_agents} agents, {dur} yr)')
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    outfile = _output_path('age_gap_distribution.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Histogram saved to {outfile}")

    print(f"\n  age_gap | mean(frac) | std(frac) |  5% CI   |  95% CI  | mean(n_rel) | std(n_rel)")
    print(f"  --------+------------+-----------+----------+----------+-------------+-----------")
    for label, mf, sf, lo, hi, mc, sd in zip(bin_labels, mean_fractions, std_fractions, ci_low, ci_high, mean_counts, std_counts):
        print(f"  {label:>7s} | {mf:10.4f} | {sf:9.4f} | {lo:8.4f} | {hi:8.4f} | {mc:11.2f} | {sd:10.2f}")

    # Overall age-gap summary (pooled across all runs, plus across-run aggregates)
    pooled = np.concatenate(all_gaps) if all_gaps else np.array([], dtype=float)
    pooled_mean = float(pooled.mean()) if pooled.size else float('nan')
    pooled_std = float(pooled.std()) if pooled.size else float('nan')
    mean_of_run_means = float(np.nanmean(per_run_gap_mean))
    std_of_run_means = float(np.nanstd(per_run_gap_mean))
    mean_of_run_stds = float(np.nanmean(per_run_gap_std))
    std_of_run_stds = float(np.nanstd(per_run_gap_std))
    print(f"\n  Age-gap summary (male - female, years):")
    print(f"    pooled across all runs : mean = {pooled_mean:7.3f}   std = {pooled_std:7.3f}   "
          f"n = {pooled.size}")
    print(f"    across-run aggregates  : mean(per-run mean) = {mean_of_run_means:7.3f} "
          f"(std {std_of_run_means:5.3f})   mean(per-run std) = {mean_of_run_stds:7.3f} "
          f"(std {std_of_run_stds:5.3f})")

    return {'mean_fractions': mean_fractions, 'std_fractions': std_fractions,
            'ci_low': ci_low, 'ci_high': ci_high,
            'mean_counts': mean_counts, 'std_counts': std_counts,
            'bin_labels': bin_labels, 'bin_edges': bin_edges,
            'per_run_counts': per_run_counts, 'per_run_fractions': per_run_fractions,
            'per_run_gap_mean': per_run_gap_mean, 'per_run_gap_std': per_run_gap_std,
            'pooled_mean': pooled_mean, 'pooled_std': pooled_std}


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_age_differences()
    test_algorithm_comparison()
    test_n_partners_distribution()
    test_n_partners_pyramid()
    test_age_gap_distribution()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
