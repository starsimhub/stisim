"""
Tests for stisim.PartnershipFormationAnalyzer.

Runs a small STI simulation, prints tabular partnership-formation counts per
gender / age bin / timestep, and produces two plots:
  1. Total male vs female new partnerships per timestep.
  2. Stacked per-age-bin subplots (youngest at bottom, oldest at top), each
     with male/female lines.
"""

import os
import shutil

import numpy as np
import pandas as pd
import pylab as pl
import starsim as ss
import stisim as sti

# Directory holding this test file. Plots are saved inside a dedicated
# subdirectory below it, which is recreated fresh on each run.
TEST_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(TEST_DIR, 'partnership_formation_analyzer')


# Default test configuration (configurable via function arguments)
DEFAULT_N_AGENTS  = 10000
DEFAULT_DUR_YEARS = 5
DEFAULT_N_SEEDS   = 3
DEFAULT_AGE_BINS  = None  # None -> PartnershipFormationAnalyzer's own default

# Plot colors per the instructions
MALE_COLOR   = 'blue'
FEMALE_COLOR = 'red'


ANALYZER_NAME = 'partnershipformationanalyzer'


def make_sim(n_agents=DEFAULT_N_AGENTS, dur=DEFAULT_DUR_YEARS, age_bins=None, rand_seed=1):
    """Construct a sim with a StructuredSexual network and the analyzer.

    Returns the sim only — fetch the analyzer post-init via
    ``sim.analyzers[ANALYZER_NAME]`` because sti.Sim deep-copies the analyzer
    during init, leaving any pre-init reference stale.
    """
    sim = sti.Sim(
        n_agents=n_agents,
        dur=dur,
        rand_seed=rand_seed,
        diseases=[sti.HIV(init_prev=0.05)],
        networks=[sti.StructuredSexual()],
        analyzers=[sti.PartnershipFormationAnalyzer(age_bins=age_bins)],
    )
    return sim


def build_dataframe(sim, analyzer):
    """Return a long-format DataFrame of (ti, year, sex, age_bin, count).

    Sums per-network counts so the resulting DataFrame is an aggregate across
    every tracked network.
    """
    yearvec = sim.t.yearvec
    rows = []
    for ti, year in enumerate(yearvec):
        for bin_str in analyzer.age_bins:
            for sex in ('f', 'm'):
                per_network = analyzer.results[f'new_partnerships_{sex}']
                cnt = int(sum(per_network[nw_safe][bin_str][ti]
                              for nw_safe in per_network))
                rows.append({'ti': ti, 'year': float(year), 'sex': sex,
                             'age_bin': bin_str, 'count': cnt})
    return pd.DataFrame(rows)


def print_tables(df, analyzer):
    """Print one wide table per sex: rows = age bins, columns = timesteps.

    Cell values are means across seeds, formatted to 2 decimal places.
    """
    for sex in ('f', 'm'):
        label = 'FEMALE' if sex == 'f' else 'MALE'
        sub = df[df.sex == sex]
        wide = sub.pivot(index='age_bin', columns='ti', values='count')
        wide = wide.reindex(analyzer.age_bins)  # keep age order
        print(f'\n=== Mean new partnerships per timestep ({label}) ===')
        with pd.option_context('display.max_columns', None,
                               'display.width', 200,
                               'display.max_rows', None,
                               'display.float_format', '{:.2f}'.format):
            print(wide)


def _title_suffix(dur, n_agents, n_seeds):
    return f'dur={dur}y, n_agents={n_agents}, n_seeds={n_seeds}'


def plot_totals(df, sim, dur, n_agents, n_seeds):
    """Total new partnerships per timestep (mean across seeds), one line per sex."""
    yearvec = sim.t.yearvec
    n_ti = len(yearvec)
    male   = (df[df.sex == 'm'].groupby('ti')['count'].sum().reindex(range(n_ti)).fillna(0).values)
    female = (df[df.sex == 'f'].groupby('ti')['count'].sum().reindex(range(n_ti)).fillna(0).values)

    # Integer year marks that fall within the simulated yearvec range; these
    # become thin black vertical gridlines on every subplot.
    year_marks = list(range(int(np.floor(yearvec[0])), int(np.floor(yearvec[-1])) + 1))

    fig, ax = pl.subplots(figsize=(9, 5))
    ax.plot(yearvec, male,   color=MALE_COLOR,   label='Male')
    ax.plot(yearvec, female, color=FEMALE_COLOR, label='Female')
    for y in year_marks:
        ax.axvline(y, color='black', linewidth=0.5)
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean number of partnerships formed')
    ax.set_title(f'Total partnerships formed per timestep\n'
                 f'({_title_suffix(dur, n_agents, n_seeds)})')
    ax.legend()
    pl.tight_layout()
    return fig


def plot_by_age_bin(df, analyzer, sim, dur, n_agents, n_seeds):
    """Stacked subplots (one per age bin), youngest at bottom, oldest at top.

    Each subplot shows the mean across seeds for male/female partner-side counts.
    """
    yearvec = sim.t.yearvec
    n_bins = len(analyzer.age_bins)
    fig, axes = pl.subplots(n_bins, 1, figsize=(9, max(1.2 * n_bins, 4)),
                            sharex=True, sharey=True)
    if n_bins == 1:
        axes = [axes]

    # Integer year marks that fall within the simulated yearvec range; these
    # become thin black vertical gridlines on every subplot.
    year_marks = list(range(int(np.floor(yearvec[0])), int(np.floor(yearvec[-1])) + 1))

    for i, bin_str in enumerate(analyzer.age_bins):
        # axes[0] is the TOP row. Youngest -> bottom -> last axis.
        ax = axes[n_bins - 1 - i]
        sub = df[df.age_bin == bin_str].sort_values('ti')
        male   = sub[sub.sex == 'm']['count'].values
        female = sub[sub.sex == 'f']['count'].values
        ax.plot(yearvec, male,   color=MALE_COLOR,   label='Male')
        ax.plot(yearvec, female, color=FEMALE_COLOR, label='Female')
        for y in year_marks:
            ax.axvline(y, color='black', linewidth=0.5)
        ax.set_ylabel(bin_str, fontsize=9, rotation=0, ha='right', va='center')

    axes[-1].set_xlabel('Year')
    axes[0].set_title(f'Mean partnerships formed per timestep by age bin '
                      f'(youngest at bottom)\n'
                      f'({_title_suffix(dur, n_agents, n_seeds)})')
    axes[0].legend(loc='upper right', fontsize=8)
    pl.tight_layout()
    return fig


def test_partnership_formation_analyzer(n_agents=DEFAULT_N_AGENTS,
                                        dur=DEFAULT_DUR_YEARS,
                                        age_bins=DEFAULT_AGE_BINS,
                                        n_seeds=DEFAULT_N_SEEDS,
                                        show_plots=False):
    """End-to-end test: run n_seeds sims in parallel, average their
    per-(sex, age_bin, ti) counts, then print tables and produce plots."""
    # Build one sim per seed and run them in parallel via ss.parallel (a thin
    # wrapper around ss.MultiSim). The returned MultiSim holds the run sims.
    sims = [make_sim(n_agents=n_agents, dur=dur, age_bins=age_bins, rand_seed=seed)
            for seed in range(1, n_seeds + 1)]
    msim = ss.parallel(*sims)
    run_sims = list(msim.sims)

    # Build per-seed long-format DataFrames and concat with a seed index, then
    # average across seeds per (ti, year, sex, age_bin) cell.
    dfs = []
    for seed_i, sim in enumerate(run_sims, start=1):
        analyzer = sim.analyzers['partnershipformationanalyzer']
        df_seed = build_dataframe(sim, analyzer)
        df_seed['seed'] = seed_i
        dfs.append(df_seed)
    df_all = pd.concat(dfs, ignore_index=True)
    df_mean = (df_all
               .groupby(['ti', 'year', 'sex', 'age_bin'], as_index=False)['count']
               .mean())

    # Use the last run sim/analyzer for shared axis info (identical structure
    # across seeds: same yearvec, same age_bins).
    sim_ref = run_sims[-1]
    analyzer_ref = sim_ref.analyzers['partnershipformationanalyzer']

    # ----- Sanity checks (on the mean DataFrame) -----
    total_mean = float(df_mean['count'].sum())
    assert total_mean > 0, 'Expected at least some new partnerships in the mean run.'

    # StructuredSexual is male-female only, so male and female means should
    # match. Allow a small tolerance for partners aging above the top bin.
    male_total_mean   = float(df_mean[df_mean.sex == 'm']['count'].sum())
    female_total_mean = float(df_mean[df_mean.sex == 'f']['count'].sum())
    denom = max(male_total_mean, female_total_mean, 1.0)
    diff_share = abs(male_total_mean - female_total_mean) / denom
    assert diff_share < 0.10, (
        f'Male/female mean totals diverge by more than 10% '
        f'(m={male_total_mean:.2f}, f={female_total_mean:.2f}).'
    )

    # No counts in pre-debut age bins (0-5, 5-10) — debut age is well above 10.
    pre_debut = float(df_mean[df_mean.age_bin.isin(['0-5', '5-10'])]['count'].sum())
    assert pre_debut == 0.0, (
        f'Found {pre_debut} mean partnerships in pre-debut age bins (0-5, 5-10).'
    )

    # Later timesteps should accumulate counts (i.e. not all activity at ti=0).
    per_ti_total = df_mean.groupby('ti')['count'].sum().values
    assert per_ti_total[1:].sum() > 0, 'No new partnerships after the first timestep.'

    # ----- Reporting -----
    print(f'\nSim: n_agents={n_agents}, dur={dur}y, '
          f'{len(sim_ref.t.yearvec)} timesteps, n_seeds={n_seeds}')
    print(f'Mean total new partnerships across seeds: {total_mean:.2f}')
    print(f'  Male mean contribution:   {male_total_mean:.2f}')
    print(f'  Female mean contribution: {female_total_mean:.2f}')

    print_tables(df_mean, analyzer_ref)
    fig_totals = plot_totals(df_mean, sim_ref, dur, n_agents, n_seeds)
    fig_by_bin = plot_by_age_bin(df_mean, analyzer_ref, sim_ref, dur, n_agents, n_seeds)

    # Recreate the output subdirectory fresh, then save plots with names
    # encoding the run's configuration.
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    stem = f'partnership_formation_dur{dur}y_n{n_agents}_seeds{n_seeds}'
    totals_path = os.path.join(OUTPUT_DIR, f'{stem}_totals.png')
    by_bin_path = os.path.join(OUTPUT_DIR, f'{stem}_by_age_bin.png')
    fig_totals.savefig(totals_path, dpi=100, bbox_inches='tight')
    fig_by_bin.savefig(by_bin_path, dpi=100, bbox_inches='tight')
    print(f'Saved plots to {OUTPUT_DIR}:\n  {os.path.basename(totals_path)}\n  {os.path.basename(by_bin_path)}')

    if show_plots:
        pl.show()
    else:
        pl.close(fig_totals)
        pl.close(fig_by_bin)

    return run_sims, analyzer_ref, df_mean


def _first_bin_for_age(analyzer, age):
    """Return the index of the first age bin containing ``age``, or -1 if none."""
    for i, (lo, hi) in enumerate(analyzer.age_ranges):
        if lo <= age < hi:
            return i
    return -1


def test_count_by_age_bin_duplicates(n_agents=200):
    """Regression: each appearance of a UID in the input array contributes one
    count. This pins the duplicate-preserving behavior of ``ss.uids`` and
    ``FloatArr`` / ``BoolArr`` indexing inside ``_count_by_age_bin``.
    """
    sim = sti.Sim(
        n_agents=n_agents, dur=1,
        diseases=[sti.HIV(init_prev=0.05)],
        networks=[sti.StructuredSexual()],
        analyzers=[sti.PartnershipFormationAnalyzer()],
    )
    sim.init()
    analyzer = sim.analyzers['partnershipformationanalyzer']
    ppl = sim.people

    # Find ANY one female and ANY one male UID whose age falls inside A bin
    female_uid = None
    male_uid = None
    for uid in range(n_agents):
        if _first_bin_for_age(analyzer, float(ppl.age[uid])) < 0:
            continue
        if bool(ppl.female[uid]) and female_uid is None:
            female_uid = uid
        elif (not bool(ppl.female[uid])) and male_uid is None:
            male_uid = uid
        if female_uid is not None and male_uid is not None:
            break
    assert female_uid is not None and male_uid is not None, ('Could not find both a female and male UID in a specified age_bin')

    # create a contrived case where this arbitrary female_uid appears 3 times, and the arbitrary male_uid appears 1 time
    uid_array = np.array([female_uid, male_uid, female_uid, female_uid])
    # determine a bin containing them (individually)
    f_bin_index = _first_bin_for_age(analyzer, float(ppl.age[female_uid]))
    m_bin_index = _first_bin_for_age(analyzer, float(ppl.age[male_uid]))

    # verify the internal _count_by_age_bin function properly counts the provided contrived data (3xfemale, 1xmale)
    f_counts = analyzer._count_by_age_bin(uid_array, female=True)
    m_counts = analyzer._count_by_age_bin(uid_array, female=False)

    assert f_counts[f_bin_index] == 3, (
        f'Female UID appeared 3 times but matching bin at index {f_bin_index} count is {f_counts[f_bin_index]} '
        f'(full counts: {f_counts.tolist()}).'
    )
    assert m_counts[m_bin_index] == 1, (
        f'Male UID appeared 1 time but matching bin at index {m_bin_index} count is {m_counts[m_bin_index]} '
        f'(full counts: {m_counts.tolist()}).'
    )
    # No spurious counts elsewhere; we ONLY passed in 3xFemale + 1xMale (default bins are non-overlapping, so each agent
    # is at most single-counted)
    assert int(f_counts.sum()) == 3, f'Stray female counts exist, full counts: {f_counts.tolist()}'
    assert int(m_counts.sum()) == 1, f'Stray male counts exist, full counts: {m_counts.tolist()}'

    # Empty input → all-zero counts
    empty = np.empty(0, dtype=int)  # there is no data in this array
    assert int(analyzer._count_by_age_bin(empty, female=True).sum()) == 0
    assert int(analyzer._count_by_age_bin(empty, female=False).sum()) == 0

    print(f'test_count_by_age_bin_duplicates: female bin {f_bin_index} = 3, '
          f'male bin {m_bin_index} = 1, empty input = 0. PASSED.')
    return


def test_count_by_age_bin_overlapping(n_agents=200):
    """Regression: with overlapping age bins, each agent is counted in every
    bin whose range contains their age (not just one).
    """
    overlapping_bins = ['15-50', '40-50']  # [40, 50) falls in BOTH
    sim = sti.Sim(
        n_agents=n_agents, dur=1,
        diseases=[sti.HIV(init_prev=0.05)],
        networks=[sti.StructuredSexual()],
        analyzers=[sti.PartnershipFormationAnalyzer(age_bins=overlapping_bins)],
    )
    sim.init()
    analyzer = sim.analyzers['partnershipformationanalyzer']
    ppl = sim.people

    # sanity check to ensure bins parsed correctly
    assert analyzer.age_bins == overlapping_bins
    assert analyzer.age_ranges == [(15.0, 50.0), (40.0, 50.0)]

    # Find a female agent aged 40-49 (in BOTH bins)
    both_bins_female = None
    for uid in range(n_agents):
        age = float(ppl.age[uid])
        if bool(ppl.female[uid]) and 40 <= age < 50:
            both_bins_female = uid
            break
    assert both_bins_female is not None, 'Could not find a female aged 40-49 as needed by this test'

    # Single occurrence: expect count of 1 in BOTH bins
    counts = analyzer._count_by_age_bin(np.array([both_bins_female]), female=True)
    assert counts.tolist() == [1, 1], (
        f'Overlapping-bins 15-50, 40-50: agent aged 40-49 should count in both bins, '
        f'got {counts.tolist()} (expected [1, 1]).'
    )

    # Repeated occurrences should still accumulate independently in each bin
    counts_repeated = analyzer._count_by_age_bin(
        np.array([both_bins_female, both_bins_female, both_bins_female]), female=True)
    assert counts_repeated.tolist() == [3, 3], (
        f'Overlapping-bins 15-50, 40-50: agent aged 40-49 should count in both bins three times in this test, '
        f'got {counts.tolist()} (expected [3, 3]).'
    )

    # Find a female aged 15-39 who should only land in 15-50
    first_bin_only_female = None
    for uid in range(n_agents):
        age = float(ppl.age[uid])
        if bool(ppl.female[uid]) and 15 <= age < 40:
            first_bin_only_female = uid
            break
    assert first_bin_only_female is not None, 'Could not find a female aged 15-40 as needed by this test'

    counts_first_bin_only = analyzer._count_by_age_bin(np.array([first_bin_only_female]), female=True)
    assert counts_first_bin_only.tolist() == [1, 0], (
        f'Female aged 15-39 should land only in 15-50, got {counts_first_bin_only.tolist()}.'
    )
    print(f'test_count_by_age_bin_overlapping: agent aged 40-49 counted in both bins, '
          f'repeated occurrences accumulate. PASSED.')
    return


if __name__ == '__main__':
    test_count_by_age_bin_duplicates()
    test_count_by_age_bin_overlapping()
    test_partnership_formation_analyzer(show_plots=True)