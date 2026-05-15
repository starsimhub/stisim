"""
Network dynamics validation: verify that network parameters (concurrency, pair
formation, relationship duration, debut age, MSM) affect behaviour as expected.
"""
import os
import sys
import time
import warnings
import stisim as sti
import numpy as np
import sciris as sc

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

    target_age_gap = 10

    sdA = 3
    sdB = 3 # 2
    sdC = 3 # 1
    age_diff_pars_older = {'teens': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdC)],
                           'young': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdB)],
                           'adult': [(target_age_gap, sdA), (target_age_gap, sdA), (target_age_gap, sdB)]}
    network = sti.MFNetwork(age_diff_pars=age_diff_pars_older)
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

    target_age_gaps = [8]
    seeds = [1, 2, 3, 4, 5] #, 4, 5, 6, 7]# , 2, 3]# , 3, 4, 5]
    dur = 10

    algorithms = ['match_pairs_existing', 'match_pairs_two_pointer_linear', 'match_pairs_two_pointer_linear_closest'] #, 'match_pairs_target_age_multipass', ]
    algorithms = ['match_pairs_two_pointer_linear', 'match_pairs_two_pointer_linear_closest', 'match_pairs_two_pointer_linear_closest_age_bounding']
    algorithms = ['match_pairs_two_pointer_linear_closest_age_bounding', 'match_pairs_two_pointer_linear_closest_age_bounding_binary_search']

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


if __name__ == '__main__':
    test_age_differences()
    # test_age_differences_sweep()
    test_algorithm_comparison()
