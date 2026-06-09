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
import warnings

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc
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
        networks=[sti.MFNetwork()],
        analyzers=[sti.PartnershipFormationAnalyzer(age_bins=age_bins)],
    )
    return sim


def build_dataframe(sim, analyzer):
    """Return a long-format DataFrame of (ti, year, sex, age_bin, count).

    Sums the per-network per-timestep formation series so the resulting
    DataFrame is an aggregate across every tracked network.
    """
    yearvec = sim.t.yearvec
    n_ti = len(yearvec)

    # {sex: {age_bin: array(n_ti)}} summed across tracked networks.
    totals = {}
    for sex in ['f', 'm']:
        # {nw: {age_bin: array(n_ti)}}, full per-timestep series (window_months=None)
        per_network = analyzer.get_n_partnerships_formed(female=(sex == 'f'))
        agg = {bin_str: np.zeros(n_ti, dtype=int) for bin_str in analyzer.age_bins}
        for bins_by_nw in per_network.values():
            for bin_str, series in bins_by_nw.items():
                agg[bin_str] += series
        totals[sex] = agg

    rows = []
    for ti, year in enumerate(yearvec):
        for bin_str in analyzer.age_bins:
            for sex in ['f', 'm']:
                rows.append({'ti': ti, 'year': float(year), 'sex': sex,
                             'age_bin': bin_str, 'count': int(totals[sex][bin_str][ti])})
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


# ===========================================================================
# Pure-function unit tests
# ---------------------------------------------------------------------------
# Exercise the static age-binning helper directly, with no analyzer instance
# and no simulation.
# ===========================================================================

def test_bin_counts_duplicates():
    """``_bin_counts`` counts each occurrence (duplicates preserved); ages
    outside every bin are skipped; empty input -> zeros."""
    ranges = [ss.parse_age_range(b) for b in ['0-5', '5-10', '10-15']]
    ages = np.array([12.0, 7.0, 12.0, 12.0, 99.0])  # 7 -> [5,10); 12 x3 -> [10,15); 99 out
    counts = sti.PartnershipFormationAnalyzer._bin_counts(ages, ranges)
    assert counts.tolist() == [0, 1, 3], counts.tolist()
    assert sti.PartnershipFormationAnalyzer._bin_counts(np.empty(0), ranges).tolist() == [0, 0, 0]
    print('test_bin_counts_duplicates: PASSED.')
    return


def test_bin_counts_overlapping():
    """``_bin_counts``: overlapping bins each count an age in their range."""
    ranges = [ss.parse_age_range(b) for b in ['15-50', '40-50']]  # [40,50) in BOTH
    assert sti.PartnershipFormationAnalyzer._bin_counts(np.array([42.0]), ranges).tolist() == [1, 1]
    assert sti.PartnershipFormationAnalyzer._bin_counts(np.array([42.0, 42.0, 42.0]), ranges).tolist() == [3, 3]
    assert sti.PartnershipFormationAnalyzer._bin_counts(np.array([20.0]), ranges).tolist() == [1, 0]
    print('test_bin_counts_overlapping: PASSED.')
    return


# ===========================================================================
# Injection-based getter tests
# ---------------------------------------------------------------------------
# The getters (get_n_partnerships_formed, get_n_partnerships_formed_per_agent,
# get_unique_partners_active_per_agent) are pure functions of the recorded
# state: self.results['edges'], self._final_uids and self._final_ti. The helper
# below injects that state directly, so each test below can assert EXACT getter
# outputs for a hand-built edge history -- which a real (stochastic) sim cannot
# give you (it only supports loose aggregate checks). This is where uniqueness
# vs non-uniqueness, duplicate handling, the trailing-window boundary,
# active-vs-formed semantics, per-network separation, and same-sex handling are
# verified to the exact value.
#
# CAVEAT: injection bypasses step() and init_results(), so these tests do NOT
# verify that the recorded state is populated correctly during a run, nor the
# `n_formed` fast path in get_n_partnerships_formed. That wiring is covered by
# the sim-level test (test_partnership_formation_analyzer) and the
# network-selection tests (test_network_*), which call sim.init()/run(). The
# 8-tuple record format documented below is the contract between step() and
# these tests; if step() changes its record layout, update the helper to match.
# ===========================================================================

def _inject_analyzer(records, uids_by_sex, final_ti, age_bins=None, nw_name='mfnetwork'):
    """Build a PartnershipFormationAnalyzer with injected state (no sim run).

    ``records`` is either a list of per-timestep 8-tuples for a single network
    (stored under ``nw_name``), or a dict ``{network_name: [records...]}`` for
    multiple tracked networks. Each record is::

        (ti, p1, p2, formation_ti, age_p1, age_p2, female_p1, female_p2)

    Sex is explicit per endpoint, so same-sex edges are expressible. This drives
    the getters directly off ``self.results['edges']`` and the cached subject
    universes, bypassing the simulation loop.
    """
    records_by_nw = records if isinstance(records, dict) else {nw_name: records}
    ana = sti.PartnershipFormationAnalyzer(age_bins=age_bins)
    ana._tracked_names = list(records_by_nw.keys())
    ana.results = {'edges': {
        nw: [(ti, np.asarray(p1, np.int64), np.asarray(p2, np.int64),
              np.asarray(fti, np.int64), np.asarray(a1, float), np.asarray(a2, float),
              np.asarray(f1, bool), np.asarray(f2, bool))
             for (ti, p1, p2, fti, a1, a2, f1, f2) in recs]
        for nw, recs in records_by_nw.items()
    }}
    ana._final_uids = {'f': np.asarray(uids_by_sex.get('f', []), np.int64),
                       'm': np.asarray(uids_by_sex.get('m', []), np.int64)}
    ana._final_ti = final_ti
    return ana


# Male-female scenario. window_months=3, final_ti=5 -> threshold=3 (ti in {3,4,5}).
# Males 10,11 (female flag False); Females 20,21,22 (female flag True).
_MF_RECORDS = [
    # (ti, p1, p2, formation_ti, age_p1, age_p2, female_p1, female_p2)
    (2, [10],             [20],             [2],          [30],             [22],
        [False],          [True]),
    # ti=3: m10-f20 forms TWICE (duplicate pair), m10-f21 forms; m11-f22 formed earlier (fti=1).
    (3, [10, 10, 10, 11], [20, 20, 21, 22], [3, 3, 3, 1], [31, 31, 31, 41], [22, 22, 18, 27],
        [False, False, False, False], [True, True, True, True]),
    (4, [10, 11],         [20, 22],         [3, 1],       [32, 42],         [23, 28],
        [False, False],   [True, True]),
    (5, [10],             [22],             [5],          [33],             [29],
        [False],          [True]),
]
_MF_UIDS = {'m': [10, 11], 'f': [20, 21, 22]}


@sc.timer()
def test_get_unique_partners_active_per_agent():
    """Unique partners per subject over the trailing window, both sexes, with
    duplicate pairings collapsed to a single unique partner."""
    sc.heading("get_unique_partners_active_per_agent: unique partners per male/female over a window")
    ana = _inject_analyzer(_MF_RECORDS, _MF_UIDS, final_ti=5)
    male = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['mfnetwork']
    female = ana.get_unique_partners_active_per_agent(female=True, window_months=3)['mfnetwork']
    assert male.tolist() == [3, 1], male.tolist()        # m10->{20,21,22}, m11->{22}
    assert female.tolist() == [1, 1, 2], female.tolist()  # f20->{10}, f21->{10}, f22->{11,10}

    # window_months=None -> whole run. Here it matches the windowed result, since
    # the only pre-window edge (ti=2 m10-f20) adds no NEW unique partner.
    male_all = ana.get_unique_partners_active_per_agent(female=False)['mfnetwork']
    female_all = ana.get_unique_partners_active_per_agent(female=True)['mfnetwork']
    assert male_all.tolist() == [3, 1], male_all.tolist()
    assert female_all.tolist() == [1, 1, 2], female_all.tolist()
    print(f'test_get_unique_partners_active_per_agent: male={male.tolist()}, female={female.tolist()}; '
          f'whole-run male={male_all.tolist()}. PASSED.')
    return


@sc.timer()
def test_get_n_partnerships_formed_per_agent():
    """Non-unique formations per agent over the window; duplicate pairings each
    count; formations outside the window are excluded."""
    sc.heading("get_n_partnerships_formed_per_agent: non-unique formations per agent over a window")
    ana = _inject_analyzer(_MF_RECORDS, _MF_UIDS, final_ti=5)
    male = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['mfnetwork']
    female = ana.get_n_partnerships_formed_per_agent(female=True, window_months=3)['mfnetwork']
    assert male.tolist() == [4, 0], male.tolist()         # m10: 3 @ ti3 + 1 @ ti5; m11: 0 (formed ti1)
    assert female.tolist() == [2, 1, 1], female.tolist()   # f20:2, f21:1, f22:1

    # window_months=None -> whole run, which additionally counts the ti=2 formation
    # (m10-f20) excluded by the window above.
    male_all = ana.get_n_partnerships_formed_per_agent(female=False)['mfnetwork']
    female_all = ana.get_n_partnerships_formed_per_agent(female=True)['mfnetwork']
    assert male_all.tolist() == [5, 0], male_all.tolist()       # m10: ti2(1)+ti3(3)+ti5(1)
    assert female_all.tolist() == [3, 1, 1], female_all.tolist()  # f20: ti2(1)+ti3(2)
    print(f'test_get_n_partnerships_formed_per_agent: male={male.tolist()}, female={female.tolist()}; '
          f'whole-run male={male_all.tolist()}, female={female_all.tolist()}. PASSED.')
    return


@sc.timer()
def test_get_n_partnerships_formed_binned():
    """Binned by age at formation. Returns {nw: {bin: array}} with the inner
    array a per-ti series (window_months=None) or a single window sum. In-window
    formation ages: male [31,31,31,33], female [22,22,18,29]."""
    sc.heading("get_n_partnerships_formed (binned): per-timestep and windowed, overlapping bins")
    ana = _inject_analyzer(_MF_RECORDS, _MF_UIDS, final_ti=5)

    # --- windowed: single trailing-window sum per bin (length-1 inner arrays) ---
    disjoint = ['15-30', '30-45']
    male = ana.get_n_partnerships_formed(female=False, age_bins=disjoint, window_months=3)['mfnetwork']
    female = ana.get_n_partnerships_formed(female=True, age_bins=disjoint, window_months=3)['mfnetwork']
    assert male['15-30'].tolist() == [0] and male['30-45'].tolist() == [4], male
    assert female['15-30'].tolist() == [4] and female['30-45'].tolist() == [0], female
    # Disjoint bins partition the observed ages -> window sums equal the per-agent totals.
    male_pa = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['mfnetwork']
    assert male['15-30'][0] + male['30-45'][0] == male_pa.sum()

    # overlapping bins each count
    overlapping = ['15-45', '30-45']  # 15-45 contains 30-45
    male_o = ana.get_n_partnerships_formed(female=False, age_bins=overlapping, window_months=3)['mfnetwork']
    female_o = ana.get_n_partnerships_formed(female=True, age_bins=overlapping, window_months=3)['mfnetwork']
    assert male_o['15-45'].tolist() == [4] and male_o['30-45'].tolist() == [4], male_o
    assert female_o['15-45'].tolist() == [4] and female_o['30-45'].tolist() == [0], female_o

    # --- per-timestep (window_months=None): inner arrays indexed by ti (len n_ti=6) ---
    # The full-run series includes ti=2 (outside the window above): male age 30 (30-45).
    male_ts = ana.get_n_partnerships_formed(female=False, age_bins=disjoint)['mfnetwork']
    # ti2 -> 1 male age 30; ti3 -> 3 males age 31; ti5 -> 1 male age 33 (all 30-45); none in 15-30.
    assert male_ts['30-45'].tolist() == [0, 0, 1, 3, 0, 1], male_ts['30-45'].tolist()
    assert male_ts['15-30'].tolist() == [0, 0, 0, 0, 0, 0], male_ts['15-30'].tolist()
    print(f'test_get_n_partnerships_formed_binned: windowed m={ {k: v.tolist() for k, v in male.items()} }; '
          f'per-ti 30-45={male_ts["30-45"].tolist()}. PASSED.')
    return


@sc.timer()
def test_same_sex_network():
    """Same-sex (male-male) edges: both endpoints are classified male via the
    stored per-endpoint sex, so unique and formed counts work with no p1/p2
    sex assumption."""
    sc.heading("Same-sex network: male-male edges counted for the male subject")
    # Males 10,11,12; all edges male-male (a male-male network, e.g. MSM).
    # final_ti=5, window=3 -> ti in {3,4,5}.
    records = [
        (3, [10, 10], [11, 12], [3, 3], [31, 31], [29, 27], [False, False], [False, False]),
        (4, [10],     [11],     [3],    [32],     [30],     [False],        [False]),  # persists
    ]
    ana = _inject_analyzer(records, {'m': [10, 11, 12], 'f': []}, final_ti=5, nw_name='msmnet')

    uniq = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['msmnet']
    assert uniq.tolist() == [2, 1, 1], uniq.tolist()      # 10->{11,12}, 11->{10}, 12->{10}
    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['msmnet']
    assert formed.tolist() == [2, 1, 1], formed.tolist()  # 10 in 2 edges @ ti3; 11,12 in 1 each
    # No females -> empty subject universe.
    assert len(ana.get_unique_partners_active_per_agent(female=True, window_months=3)['msmnet']) == 0
    print(f'test_same_sex_network: unique={uniq.tolist()}, formed={formed.tolist()}. PASSED.')
    return


@sc.timer()
def test_active_vs_formed():
    """Active-vs-formed distinction: a partnership formed BEFORE the window but
    still active during it is counted by get_unique_partners_active_per_agent
    (active edges) but NOT by get_n_partnerships_formed_per_agent (formed in
    window)."""
    sc.heading("Active-vs-formed: pre-window partnership counts as active, not as formed")
    # final_ti=5, window=3 -> threshold=3 (ti in {3,4,5}).
    # m10-f20 formed at ti=1 (before window) but active throughout it.
    # m11-f21 formed at ti=4 (inside the window).
    records = [
        (3, [10],     [20],     [1],    [31],     [22],     [False],        [True]),
        (4, [10, 11], [20, 21], [1, 4], [32, 41], [23, 19], [False, False], [True, True]),
        (5, [10, 11], [20, 21], [1, 4], [33, 42], [24, 20], [False, False], [True, True]),
    ]
    ana = _inject_analyzer(records, {'m': [10, 11], 'f': [20, 21]}, final_ti=5)

    active = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['mfnetwork']
    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['mfnetwork']
    assert active.tolist() == [1, 1], active.tolist()   # m10->{20} active; m11->{21} active
    assert formed.tolist() == [0, 1], formed.tolist()   # m10 formed pre-window=0; m11 formed @ ti4=1
    print(f'test_active_vs_formed: active={active.tolist()} vs formed={formed.tolist()}. PASSED.')
    return


@sc.timer()
def test_multiple_networks():
    """Two tracked networks return independent per-network results; a single
    agent (uid=2) present in both has independent partnerships per network."""
    sc.heading("Multiple networks: independent per-network results; shared agent uid=2")
    # TODO: consider adding a real MF+SW integration smoke test (a sti.Sim with
    # both MFNetwork and SWNetwork) to complement this injection-based check.
    # final_ti=5, window=3 -> threshold=3.
    records_by_nw = {
        'mfnetwork': [
            (3, [2], [20], [3], [31], [22], [False], [True]),    # uid2 forms with f20 in MF
        ],
        'swnetwork': [
            (3, [2], [30], [3], [31], [25], [False], [True]),    # uid2 (client) forms with fsw30
            (4, [2], [31], [4], [32], [26], [False], [True]),    # uid2 forms with fsw31
        ],
    }
    uids = {'m': [2], 'f': [20, 30, 31]}
    ana = _inject_analyzer(records_by_nw, uids, final_ti=5)

    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)
    assert set(formed.keys()) == {'mfnetwork', 'swnetwork'}, list(formed.keys())
    assert formed['mfnetwork'].tolist() == [1], formed['mfnetwork'].tolist()   # uid2: 1 in MF
    assert formed['swnetwork'].tolist() == [2], formed['swnetwork'].tolist()   # uid2: 2 in SW

    uniq = ana.get_unique_partners_active_per_agent(female=False, window_months=3)
    assert uniq['mfnetwork'].tolist() == [1], uniq['mfnetwork'].tolist()       # {20}
    assert uniq['swnetwork'].tolist() == [2], uniq['swnetwork'].tolist()       # {30,31}

    # Female side stays separated per network too.
    f_formed = ana.get_n_partnerships_formed_per_agent(female=True, window_months=3)
    assert f_formed['mfnetwork'].tolist() == [1, 0, 0], f_formed['mfnetwork'].tolist()  # f20 only
    assert f_formed['swnetwork'].tolist() == [0, 1, 1], f_formed['swnetwork'].tolist()  # f30,f31
    print(f"test_multiple_networks: uid2 formed MF={formed['mfnetwork'].tolist()}, "
          f"SW={formed['swnetwork'].tolist()}. PASSED.")
    return


# ===========================================================================
# Construction- and sim-level behavior tests (no injected state)
# ---------------------------------------------------------------------------
# Build a real analyzer (and, where needed, a real sim) to exercise __init__
# validation and init_results() network selection / warnings / errors.
# ===========================================================================

@sc.timer()
def test_duplicate_age_bins_rejected():
    """Duplicate age-bin strings collide as dict keys and must be rejected at
    construction, with every duplicated value reported."""
    sc.heading("Duplicate age_bins raise ValueError listing all duplicates")
    try:
        sti.PartnershipFormationAnalyzer(age_bins=['0-5', '5-10', '0-5', '5-10', '10-15'])
    except ValueError as e:
        msg = str(e)
        assert '0-5' in msg and '5-10' in msg, msg
        assert '10-15' not in msg, msg  # non-duplicated bin must not be listed
        print(f'test_duplicate_age_bins_rejected: raised as expected -> {msg}. PASSED.')
        return
    raise AssertionError('Expected ValueError for duplicate age_bins, none raised.')


@sc.timer()
def test_network_missing_raises():
    """Naming a network absent from the sim raises ValueError during init."""
    sc.heading("networks=[<missing>] raises ValueError")
    sim = sti.Sim(n_agents=100, dur=1, diseases=[sti.HIV(init_prev=0.05)],
                  networks=[sti.MFNetwork()],
                  analyzers=[sti.PartnershipFormationAnalyzer(networks=['bogusnet'])])
    try:
        sim.init()
    except ValueError as e:
        assert 'bogusnet' in str(e), str(e)
        print(f'test_network_missing_raises: raised as expected -> {e}. PASSED.')
        return
    raise AssertionError('Expected ValueError for a missing network, none raised.')


@sc.timer()
def test_network_nonsexual_warns_and_skips():
    """A named non-sexual network warns and is skipped; the default (networks=None)
    selects only sexual networks (the non-sexual one is excluded with no warning)."""
    sc.heading("Non-sexual network warns + skipped; default None selects only sexual nets")

    # Explicitly naming the non-sexual maternal network -> warning + skip.
    sim = sti.Sim(n_agents=100, dur=1, diseases=[sti.HIV(init_prev=0.05)],
                  networks=[sti.MFNetwork(), ss.MaternalNet()],
                  analyzers=[sti.PartnershipFormationAnalyzer(networks=['maternalnet'])])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sim.init()
    ana = sim.analyzers['partnershipformationanalyzer']
    assert any('maternalnet' in str(w.message) for w in caught), [str(w.message) for w in caught]
    assert len(ana._tracked_names) == 0, ana._tracked_names

    # Default None -> only the sexual network is tracked (no maternalnet, no warning).
    sim2 = sti.Sim(n_agents=100, dur=1, diseases=[sti.HIV(init_prev=0.05)],
                   networks=[sti.MFNetwork(), ss.MaternalNet()],
                   analyzers=[sti.PartnershipFormationAnalyzer()])
    sim2.init()
    ana2 = sim2.analyzers['partnershipformationanalyzer']
    assert ana2._tracked_names == ['mfnetwork'], ana2._tracked_names
    print(f'test_network_nonsexual_warns_and_skips: explicit skip -> {ana._tracked_names}, '
          f'default -> {ana2._tracked_names}. PASSED.')
    return


if __name__ == '__main__':
    # Pure-function unit tests
    test_bin_counts_duplicates()
    test_bin_counts_overlapping()

    # Injection-based getter tests (exact-value verification, no sim)
    test_get_unique_partners_active_per_agent()
    test_get_n_partnerships_formed_per_agent()
    test_get_n_partnerships_formed_binned()
    test_same_sex_network()
    test_active_vs_formed()
    test_multiple_networks()

    # Construction- and sim-level behavior tests
    test_duplicate_age_bins_rejected()
    test_network_missing_raises()
    test_network_nonsexual_warns_and_skips()

    # End-to-end integration (real multi-seed sim + plots/tables)
    test_partnership_formation_analyzer(show_plots=True)