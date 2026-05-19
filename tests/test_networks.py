"""
Network dynamics validation: verify that network parameters (concurrency, pair
formation, relationship duration, debut age, MSM) affect behaviour as expected.
"""
import os
import sys
import time
import warnings
import stisim as sti
import starsim as ss
import numpy as np
import sciris as sc
from collections import defaultdict

# Allow `from test_network_proposal import _match_pairs_fast` whether pytest
# collects from the project root or this file is run directly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def test_age_differences(n_agents=100000):
    """
    Sexual network with default relationship type distribution: run for 1 month and
    analyze partner age differences. Produces a scatterplot of female vs male partner
    ages with a best-fit regression line, and prints the mean and standard deviation
    of the age gap.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    target_age_gap = 5

    sdA = 3
    sdB = 3 # 2
    sdC = 3 # 1
    age_diff_pars_older = {'teens': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdC)],
                           'young': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdB)],
                           'adult': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdB)]}
    network = sti.MFNetwork(age_diff_pars=age_diff_pars_older,
                            match_pairs_algo='match_pairs_two_pointer_linear_closest_age_bounding_binary_search')
    sim = sti.Sim(n_agents=n_agents, networks=[network], dur=1, rand_seed=1)
    sim.run()

    nw = sim.networks.mfnetwork
    male_ages = np.array(nw.edges.age_p1, dtype=float)
    female_ages = np.array(nw.edges.age_p2, dtype=float)

    # Restrict both axes to 15–75 years
    min_age = 0
    max_age = 1000
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
    ax.set_title('Relationship partner ages')
    ax.set_xlim(min_age, max_age)
    ax.set_ylim(min_age, max_age)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outfile = sc.thisdir(__file__, 'relationship_age_differences.png')
    plt.savefig(outfile, dpi=100)
    plt.close(fig)
    print(f"Scatterplot saved to {outfile}")

    return sim



def _run_age_diff_sim(target_age_gap, seed, algorithm, n_agents=10000, dur=1):
    age_diff_pars = {
        'teens': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
        'young': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
        'adult': [(target_age_gap, 3), (target_age_gap, 3), (target_age_gap, 3)],
    }
    network = sti.MFNetwork(age_diff_pars=age_diff_pars, match_pairs_algo=algorithm)
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


# def test_age_differences_sweep(n_agents=10000):
#     """Sweep ``target_age_gap`` in [2,4,6,8,10,12,14,16] x seeds [1,2,3].
#
#     Runs the current :meth:`MFNetwork.match_pairs` for each combo and reports
#     the mean reported partnership age difference and mean reported standard
#     deviation per ``target_age_gap``. Useful for characterizing how the
#     matching algorithm reproduces the requested age-gap distribution.
#     """
#     target_age_gaps = [4, 6, 8, 10]
#     seeds = [1, 2, 3, 4, 5]
#     print(f"Sweep: {len(target_age_gaps)} target_gaps x {len(seeds)} seeds "
#           f"= {len(target_age_gaps)*len(seeds)} sims, n_agents={n_agents}")
#     results = _run_age_diff_sweep(target_age_gaps=target_age_gaps,
#                                   seeds=seeds,
#                                   n_agents=n_agents,
#                                   label='match_pairs')
#     _print_sweep_summary(results, label='match_pairs (current)')
#     return results



def test_algorithm_comparison(n_agents=100000):

    target_age_gaps = [5]
    seeds = [1, 2, 3, 4, 5] #, 4, 5, 6, 7]# , 2, 3]# , 3, 4, 5]
    dur = 10

    algorithms = ['match_pairs_existing', 'match_pairs_two_pointer_linear', 'match_pairs_two_pointer_linear_closest'] #, 'match_pairs_target_age_multipass', ]
    algorithms = ['match_pairs_two_pointer_linear', 'match_pairs_two_pointer_linear_closest', 'match_pairs_two_pointer_linear_closest_age_bounding']
    algorithms = ['match_pairs_two_pointer_linear_closest_age_bounding', 'match_pairs_two_pointer_linear_closest_age_bounding_binary_search']
    algorithms = ['match_pairs_two_pointer_linear_closest_age_bounding_binary_search']

    all_results = {}
    for algo_name in algorithms:
        print(f"\n=== Sweep: {algo_name} ===")
        all_results[algo_name] = _run_age_diff_sweep(target_age_gaps, seeds, algorithm=algo_name,
                                                     n_agents=n_agents, dur=dur, label=algo_name)
        _print_sweep_summary(all_results[algo_name], label=algo_name)

    # Side-by-side: mean age diff (mean over seeds), one column per algorithm
    print("\n=== Side-by-side: mean age diff per target_gap ===")
    header = "  target_gap | " + " | ".join(f"{n:>30s}" for n in algorithms)
    sep    = "  " + "-"*9 + "-+-" + "-+-".join(["-"*30]*len(algorithms))
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(np.nanmean(all_results[n][target_age_gap]['means'])):30.3f}"
                 for n in algorithms]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean realized std per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(np.nanmean(all_results[n][target_age_gap]['stds'])):30.3f}"
                 for n in algorithms]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean partnership count per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{int(np.mean(all_results[n][target_age_gap]['n_pairs'])):30d}"
                 for n in algorithms]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: mean run time (s) per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(all_results[n][target_age_gap]['time_mean']):30.3f}"
                 for n in algorithms]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    print("\n=== Side-by-side: run time std (s) per target_gap ===")
    print(header)
    print(sep)
    for target_age_gap in target_age_gaps:
        cells = [f"{float(all_results[n][target_age_gap]['time_stds']):30.3f}"
                 for n in algorithms]
        print(f"  {target_age_gap:9d} | " + " | ".join(cells))

    return {'stats': all_results}


class NPartnersAnalyzer(ss.Analyzer):
    """Track unique female partners per male over a trailing month window.

    Walks the active MFNetwork edges each step and remembers, per male, the
    most recent timestep on which each female partner appeared. After the sim
    finishes, ``get_partner_counts`` returns the number of unique female
    partners each male saw in the trailing ``window_months`` months.

    Args:
        network (str): Network name to inspect (default ``'mfnetwork'``).
        window_months (int): Trailing window in months; assumes monthly dt
            (one ti per month) so this is also the number of trailing steps.
    """

    def __init__(self, network='mfnetwork', window_months=12, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.network = network
        self.window_months = window_months
        self._records = []  # list of (ti, male_uids_array, female_uids_array)
        self._final_male_uids = None
        self._final_ti = None
        return

    def step(self):
        sim = self.sim
        nw = sim.networks[self.network]
        p1 = np.asarray(nw.edges.p1, dtype=np.int64)
        p2 = np.asarray(nw.edges.p2, dtype=np.int64)
        if len(p1):
            self._records.append((self.ti, p1, p2))
        # Cache so we can compute counts even after sim is finalized/shrunk
        self._final_male_uids = np.asarray(sim.people.male.uids, dtype=np.int64)
        self._final_ti = self.ti
        return

    def get_partner_counts(self):
        """Return per-male unique-partner counts over the trailing window."""
        threshold = self._final_ti - self.window_months + 1
        partner_sets = defaultdict(set)
        for ti, males, females in self._records:
            if ti < threshold:
                continue
            for m, f in zip(males, females):
                partner_sets[int(m)].add(int(f))

        male_uids = self._final_male_uids
        counts = np.zeros(len(male_uids), dtype=np.int64)
        for i, uid in enumerate(male_uids):
            counts[i] = len(partner_sets.get(int(uid), ()))
        return counts


def test_n_partners_distribution(n_agents=25000, n_runs=5, dur=25, window_months=12, target_age_gap=5):
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
        network = sti.MFNetwork(age_diff_pars=age_diff_pars,
                                match_pairs_algo='match_pairs_two_pointer_linear_closest_age_bounding_binary_search')
        analyzer = NPartnersAnalyzer(network='mfnetwork', window_months=window_months)
        sim = sti.Sim(n_agents=n_agents, networks=[network], analyzers=[analyzer],
                      dur=dur, rand_seed=seed)
        sim.run()
        # sti.Sim copies the analyzer during init, so retrieve the live instance
        sim_analyzer = sim.analyzers.npartnersanalyzer
        counts = sim_analyzer.get_partner_counts()
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
    outfile = sc.thisdir(__file__, 'n_partners_distribution.png')
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


def test_age_gap_distribution(n_agents=25000, n_runs=5, dur=25, window_months=12,
                              target_age_gap=5, gap_bin_width=1):
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
    max_allowed_age_delta = target_age_gap + 2 * sdA


    seeds = list(range(1, n_runs + 1))
    # Integer-centered bins from -gap_range to +gap_range inclusive
    n_bins = int(round(2 * max_allowed_age_delta / gap_bin_width)) + 1
    half = gap_bin_width / 2
    bin_edges = np.linspace(-(max_allowed_age_delta+1) - half, (max_allowed_age_delta+1) + half, n_bins + 1)
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
        network = sti.MFNetwork(age_diff_pars=age_diff_pars,
                                match_pairs_algo='match_pairs_two_pointer_linear_closest_age_bounding_binary_search')

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
    outfile = sc.thisdir(__file__, 'age_gap_distribution.png')
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
    test_age_differences()
    # test_age_differences_sweep()
    test_algorithm_comparison()
    test_n_partners_distribution()
    test_age_gap_distribution()
