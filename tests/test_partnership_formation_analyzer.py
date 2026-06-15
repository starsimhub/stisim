"""
Tests for stisim.PartnershipFormationAnalyzer.

Runs a small STI simulation, prints tabular partnership-formation counts per
gender / age bin / timestep, and produces two plots:
  1. Total male vs female new partnerships per timestep.
  2. Stacked per-age-bin subplots (youngest at bottom, oldest at top), each
     with male/female lines.
"""


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pylab as pl
import sciris as sc
import shutil
import starsim as ss
import stisim as sti
import warnings


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


def make_sim(n_agents=DEFAULT_N_AGENTS, dur=DEFAULT_DUR_YEARS, age_bins=None, rand_seed=1):
    """Construct a sim with an MFNetwork and the analyzer.

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
                                        n_seeds=DEFAULT_N_SEEDS):
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

    # MFNetwork is male-female only, so male and female means should
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

    # TODO: perhaps move these plots after some review into devtest for this relationships
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

    return run_sims, analyzer_ref, df_mean


# ===========================================================================
# Injection-based getter tests
# ---------------------------------------------------------------------------
# The getters (get_n_partnerships_formed, get_n_partnerships_formed_per_agent,
# get_unique_partners_active_per_agent) are pure functions of the recorded
# state: self.results['relationships'] (one interval row per edge), self._final_uids,
# and self._final_ti. The helper below injects that state directly, so each
# test asserts EXACT getter outputs for a hand-built edge history -- which a
# real (stochastic) sim cannot give (it only supports loose aggregate checks).
# This is where uniqueness vs non-uniqueness, the trailing-window boundary
# (half-open [ti_formed, ti_expired)), active-vs-formed semantics,
# per-network separation, and same-sex handling are verified to the value.
#
# CAVEAT: injection bypasses step() and init_results(), so these tests do NOT
# verify that edge intervals are recorded correctly during a run (formation
# detection + ti_expired patching from expired_this_ti). That wiring is covered
# by the sim-level test (test_partnership_formation_analyzer) and the
# construction/sim-level tests below (which call sim.init()/run()). The edge-row
# layout below is the contract with step(); if it changes, update the helper.
# ===========================================================================

def _inject_analyzer(edges_by_nw, uids_by_sex, final_ti, age_bins=None):
    """Build a PartnershipFormationAnalyzer with injected edge intervals (no sim).

    ``edges_by_nw`` is a dict ``{network_name: [rows...]}`` of the tracked
    networks' edge tables. Each row is::

        [p1, p2, ti_formed, ti_expired, age_p1, age_p2, female_p1, female_p2]

    ``ti_expired`` is the expiry timestep (active interval is the half-open
    ``[ti_formed, ti_expired)``), or ``None`` if the edge is still active at
    run end. Sex is explicit per endpoint, so same-sex edges are expressible.
    This drives the getters directly off ``self.results['relationships']`` and the
    cached surviving-agent uids, bypassing the simulation loop.
    """
    ana = sti.PartnershipFormationAnalyzer(age_bins=age_bins)
    ana._tracked_nw_names = list(edges_by_nw.keys())
    ana.results = {'relationships': {nw: [list(row) for row in rows] for nw, rows in edges_by_nw.items()}}
    # ana._indicies_of_rels_in_table_for_nw = {nw: {} for nw in edges_by_nw} # not strictly unless running a sim
    ana._final_uids = {'f': np.asarray(uids_by_sex.get('f', []), np.int64),
                       'm': np.asarray(uids_by_sex.get('m', []), np.int64)}
    ana._final_ti = final_ti
    return ana


# Male-female scenario as edge intervals. final_ti=5; window_months=3 ->
# threshold=3, i.e. window timesteps {3,4,5} (win_start=3, win_end=5).
# Males 10,11 (female flag False); Females 20,21,22 (female flag True).
#
# Format: [p1, p2, ti_formed, ti_expired, age_p1, age_p2, female_p1, female_p2]
# The active ti values below are known relative to usage, with "final_ti" being 5.
_MF_EDGES = [
    [10, 20, 2, 4,    30, 22, False, True],   # E1: active {2,3}
    [10, 21, 3, None, 31, 18, False, True],   # E2: active {3,4,5}
    [10, 22, 5, None, 33, 29, False, True],   # E3: active {5}
    [11, 22, 1, None, 41, 27, False, True],   # E4: formed pre-window, active throughout
    [10, 20, 5, None, 33, 24, False, True],   # E5: re-form of 10-20 (distinct edge)
    [11, 21, 0, 2,    40, 17, False, True],   # E6: active {0,1}, ends before window
]
_MF_UIDS = {'m': [10, 11], 'f': [20, 21, 22]}


@sc.timer()
def test_get_unique_partners_active_per_agent():
    """Unique partners per subject on edges ACTIVE during the window, both sexes;
    duplicate partners collapse; an edge that ended before the window is dropped."""
    sc.heading("get_unique_partners_active_per_agent: unique active partners per male/female")
    ana = _inject_analyzer({'mfnetwork': _MF_EDGES}, _MF_UIDS, final_ti=5)

    # window=3 (active during {3,4,5}). E1 (ti_expired=4 > win_start=3) is active
    # (last-active ti=3); E6 (ti_expired=2) ended before the window -> dropped.
    male = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['mfnetwork']
    female = ana.get_unique_partners_active_per_agent(female=True, window_months=3)['mfnetwork']
    assert male.tolist()   == [3, 1],    male.tolist()         # m10->{20,21,22}; m11->{22} (E6 dropped)
    assert female.tolist() == [1, 1, 2], female.tolist()  # f20->{10}; f21->{10}; f22->{10,11}

    # whole run additionally includes E6: m11->{22,21}=2 and f21->{10,11}=2.
    male_all = ana.get_unique_partners_active_per_agent(female=False)['mfnetwork']
    female_all = ana.get_unique_partners_active_per_agent(female=True)['mfnetwork']
    assert male_all.tolist()   == [3, 2],    male_all.tolist()
    assert female_all.tolist() == [1, 2, 2], female_all.tolist()
    print(f'test_get_unique_partners_active_per_agent: window male={male.tolist()}, '
          f'whole-run male={male_all.tolist()}. PASSED.')
    return


@sc.timer()
def test_get_n_partnerships_formed_per_agent():
    """Non-unique formations per agent; a pair re-forming counts each time;
    formations outside the window are excluded."""
    sc.heading("get_n_partnerships_formed_per_agent: non-unique formations per agent")
    ana = _inject_analyzer({'mfnetwork': _MF_EDGES}, _MF_UIDS, final_ti=5)

    male = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['mfnetwork']
    female = ana.get_n_partnerships_formed_per_agent(female=True, window_months=3)['mfnetwork']
    assert male.tolist()   == [3, 0],    male.tolist()         # m10: E2,E3,E5 (ft 3,5,5); m11: 0 (ft 1,0)
    assert female.tolist() == [1, 1, 1], female.tolist()  # f20:E5; f21:E2; f22:E3

    # whole run, all rels counted: m10 also gets E1(ft2) -> 4; m11 gets E4,E6 -> 2. Non-unique:
    # m10-f20 counts twice (E1+E5), so m10 formed=4 > its unique partners=3.
    male_all = ana.get_n_partnerships_formed_per_agent(female=False)['mfnetwork']
    female_all = ana.get_n_partnerships_formed_per_agent(female=True)['mfnetwork']
    assert male_all.tolist()   == [4, 2],    male_all.tolist()
    assert female_all.tolist() == [2, 2, 2], female_all.tolist()
    print(f'test_get_n_partnerships_formed_per_agent: window male={male.tolist()}, '
          f'whole-run male={male_all.tolist()}. PASSED.')
    return


@sc.timer()
def test_get_n_partnerships_formed_binned():
    """Binned by age at formation: {nw:{bin:array}}; per-ti when window_months=None,
    a single window sum otherwise. Covers disjoint, overlapping, and out-of-range bins."""
    sc.heading("get_n_partnerships_formed (binned): per-timestep, windowed, overlapping, out-of-range")
    ana = _inject_analyzer({'mfnetwork': _MF_EDGES}, _MF_UIDS, final_ti=5)

    # In-window (ft>=3) male formation ages: E2(31),E3(33),E5(33).
    # In-window female formation ages: E2 f21(18), E3 f22(29), E5 f20(24).
    disjoint = ['15-30', '30-45']
    male = ana.get_n_partnerships_formed(female=False, age_bins=disjoint, window_months=3)['mfnetwork']
    female = ana.get_n_partnerships_formed(female=True, age_bins=disjoint, window_months=3)['mfnetwork']
    assert male['15-30'].tolist()   == [0] and male['30-45'].tolist()   == [3], male
    assert female['15-30'].tolist() == [3] and female['30-45'].tolist() == [0], female

    # Overlapping bins each count the same formation (15-45 contains 30-45, which also contains all male ages).
    overlapping = ['15-45', '30-45']
    male_o = ana.get_n_partnerships_formed(female=False, age_bins=overlapping, window_months=3)['mfnetwork']
    assert male_o['15-45'].tolist() == [3] and male_o['30-45'].tolist() == [3], male_o

    # Ages outside every bin are not counted (all male ages are >= 30).
    narrow = ana.get_n_partnerships_formed(female=False, age_bins=['15-25'], window_months=3)['mfnetwork']
    assert narrow['15-25'].tolist() == [0], narrow

    # per-timestep (window_months=None): inner arrays indexed by ti (len n_ti=6).
    male_ts = ana.get_n_partnerships_formed(female=False, age_bins=disjoint)['mfnetwork']
    # 30-45 male ages by ti_formed: ft0:40(E6), ft1:41(E4), ft2:30(E1), ft3:31(E2), ft5:33,33(E3,E5).
    assert male_ts['30-45'].tolist() == [1, 1, 1, 1, 0, 2], male_ts['30-45'].tolist()
    assert male_ts['15-30'].tolist() == [0, 0, 0, 0, 0, 0], male_ts['15-30'].tolist()
    print(f"test_get_n_partnerships_formed_binned: window male={ {k: v.tolist() for k, v in male.items()} }; "
          f"per-ti 30-45={male_ts['30-45'].tolist()}. PASSED.")
    return


@sc.timer()
def test_active_vs_formed():
    """A partnership formed BEFORE the window but still active during it is counted
    by get_unique_partners_active_per_agent, but NOT by
    get_n_partnerships_formed_per_agent."""
    sc.heading("Active-vs-formed: pre-window partnership counts as active, not as formed")
    # final_ti=5, window=3 (threshold=3). m10-f20 formed ti1 (pre-window), still
    # active (ti_expired=None). m11-f21 formed ti4 (inside the window).
    edges = [
        [10, 20, 1, None, 31, 22, False, True],
        [11, 21, 4, None, 41, 19, False, True],
    ]
    ana = _inject_analyzer({'mfnetwork': edges}, {'m': [10, 11], 'f': [20, 21]}, final_ti=5)
    active = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['mfnetwork']
    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['mfnetwork']
    assert active.tolist() == [1, 1], active.tolist()   # m10 active{20}; m11 active{21}
    assert formed.tolist() == [0, 1], formed.tolist()   # m10 formed pre-window=0; m11 formed ti4=1
    print(f'test_active_vs_formed: active={active.tolist()} vs formed={formed.tolist()}. PASSED.')
    return


@sc.timer()
def test_half_open_window_boundary():
    """Half-open active interval [ti_formed, ti_expired): an edge whose
    ti_expired equals win_start is NOT active in the window (its last-active ti is
    win_start-1); a single-timestep edge inside the window is counted."""
    sc.heading("Half-open boundary: ti_expired == win_start is excluded")
    # final_ti=5, window=3 -> win_start=3, win_end=5.
    edges = [
        [10, 20, 1, 3,    31, 22, False, True],   # B1: active {1,2}; ti_expired==win_start -> excluded
        [11, 21, 1, 4,    41, 19, False, True],   # B2: active {1,2,3}; last-active 3 IS in window
        [10, 22, 4, 5,    32, 24, False, True],   # B3: active {4} (single timestep) in window
    ]
    ana = _inject_analyzer({'mfnetwork': edges}, {'m': [10, 11], 'f': [20, 21, 22]}, final_ti=5)
    win = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['mfnetwork']
    assert win.tolist() == [1, 1], win.tolist()       # m10->{22} (B1 excluded, B3 in); m11->{21}
    allrun = ana.get_unique_partners_active_per_agent(female=False)['mfnetwork']
    assert allrun.tolist() == [2, 1], allrun.tolist()  # whole run: m10->{20(B1),22(B3)}=2; m11->{21}
    print(f'test_half_open_window_boundary: window={win.tolist()}, whole-run={allrun.tolist()}. PASSED.')
    return


@sc.timer()
def test_same_sex_network():
    """Male-male edges: both endpoints are classified male via the stored
    per-endpoint sex, so unique and formed counts work with no p1/p2 assumption."""
    sc.heading("Same-sex network: male-male edges counted for the male subject")
    # Males 10,11,12; all edges male-male (e.g. an MSM network). final_ti=5, window=3.
    edges = [
        [10, 11, 3, None, 31, 29, False, False],   # m10-m11 formed ti3
        [10, 12, 3, None, 31, 27, False, False],   # m10-m12 formed ti3
    ]
    ana = _inject_analyzer({'msmnet': edges}, {'m': [10, 11, 12], 'f': []}, final_ti=5)

    uniq = ana.get_unique_partners_active_per_agent(female=False, window_months=3)['msmnet']
    assert uniq.tolist() == [2, 1, 1], uniq.tolist()      # 10->{11,12}; 11->{10}; 12->{10}
    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)['msmnet']
    assert formed.tolist() == [2, 1, 1], formed.tolist()  # both endpoints male -> 10 counts in 2 edges
    # No females -> empty surviving-agent set.
    f_formed = ana.get_unique_partners_active_per_agent(female=True, window_months=3)['msmnet']
    assert len(f_formed) == 0
    print(f'test_same_sex_network: unique={uniq.tolist()}, formed={formed.tolist()}. PASSED.')
    return


@sc.timer()
def test_multiple_networks():
    """Two tracked networks return independent per-network results; a shared
    agent (uid=2) has independent partnerships per network."""
    sc.heading("Multiple networks: independent per-network results; shared agent uid=2")
    # TODO: consider adding a real MF+SW integration smoke test (a sti.Sim with
    #  both MFNetwork and SWNetwork) to complement this injection-based check.
    #  final_ti=5, window=3, and/or StructuredSexual as a single network.
    edges_by_nw = {
        'mfnetwork': [
            [2, 20, 3, None, 31, 22, False, True],    # uid2 forms with f20 in MF
        ],
        'swnetwork': [
            [2, 30, 2, 3, 31, 25, False, True],       # uid2 forms with fsw30, outside window
            [2, 30, 3, 4, 31, 25, False, True],       # uid2 forms with fsw30
            [2, 30, 4, None, 31, 25, False, True],    # uid2 re-forms with fsw30
            [2, 31, 4, None, 32, 26, False, True],    # uid2 forms with fsw31
        ],
    }
    uids = {'m': [2], 'f': [20, 30, 31]}
    ana = _inject_analyzer(edges_by_nw, uids, final_ti=5)

    formed = ana.get_n_partnerships_formed_per_agent(female=False, window_months=3)
    assert set(formed.keys()) == {'mfnetwork', 'swnetwork'}, list(formed.keys())
    assert formed['mfnetwork'].tolist() == [1], formed['mfnetwork'].tolist()   # uid2: 1 in MF
    assert formed['swnetwork'].tolist() == [3], formed['swnetwork'].tolist()   # uid2: 2 in SW

    uniq = ana.get_unique_partners_active_per_agent(female=False, window_months=3)
    assert uniq['mfnetwork'].tolist() == [1], uniq['mfnetwork'].tolist()       # {20}
    assert uniq['swnetwork'].tolist() == [2], uniq['swnetwork'].tolist()       # {30,31}

    # Female side stays separated per network too.
    f_formed = ana.get_n_partnerships_formed_per_agent(female=True, window_months=3)
    assert f_formed['mfnetwork'].tolist() == [1, 0, 0], f_formed['mfnetwork'].tolist()  # f20 only
    assert f_formed['swnetwork'].tolist() == [0, 2, 1], f_formed['swnetwork'].tolist()  # f30,f31
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
    assert len(ana._tracked_nw_names) == 0, ana._tracked_nw_names

    # Default None -> only the sexual network is tracked (no maternalnet, no warning).
    sim2 = sti.Sim(n_agents=100, dur=1, diseases=[sti.HIV(init_prev=0.05)],
                   networks=[sti.MFNetwork(), ss.MaternalNet()],
                   analyzers=[sti.PartnershipFormationAnalyzer()])
    sim2.init()
    ana2 = sim2.analyzers['partnershipformationanalyzer']
    assert ana2._tracked_nw_names == ['mfnetwork'], ana2._tracked_nw_names
    print(f'test_network_nonsexual_warns_and_skips: explicit skip -> {ana._tracked_nw_names}, '
          f'default -> {ana2._tracked_nw_names}. PASSED.')
    return sim


@sc.timer()
def test_network_unsupported_expiry_warns_and_skips():
    """A network that cannot record all expirations (records_all_expirations=False,
    e.g. MSMScaleFreeNetwork) is skipped with a warning; other sexual networks are
    still tracked."""
    sc.heading("Network with records_all_expirations=False warns + skipped")
    sim = sti.Sim(n_agents=200, dur=1, diseases=[sti.HIV(init_prev=0.05)],
                  networks=[sti.MFNetwork(), sti.MSMScaleFreeNetwork()],
                  analyzers=[sti.PartnershipFormationAnalyzer()])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sim.init()
    ana = sim.analyzers['partnershipformationanalyzer']
    assert any('record all edge expirations' in str(w.message) for w in caught), [str(w.message) for w in caught]
    assert ana._tracked_nw_names == ['mfnetwork'], ana._tracked_nw_names
    assert 'msmscalefreenetwork' in ana._skipped_nw_names, ana._skipped_nw_names
    print(f'test_network_unsupported_expiry_warns_and_skips: tracked={ana._tracked_nw_names}, '
          f'skipped={ana._skipped_nw_names}. PASSED.')
    return sim


@sc.timer()
def test_records_expired_end_to_end(deaths=False):
    """A real run records edge expirations through the full mechanism
    (_on_edge_dissolution append + remove_uids override + step() drain +
    finalize): some edges get a finite ti_expired, finite intervals are valid
    (ti_expired > ti_formed), and the still-active (None) edges exactly match
    the network's active edges at run end."""
    sc.heading("End-to-end: MFNetwork records expirations; intervals consistent")
    if deaths:
        sim = sti.Sim(n_agents=500, dur=5, rand_seed=1, diseases=[sti.HIV(init_prev=0.05)],
                      networks=[sti.MFNetwork()], demographics=[ss.Deaths(death_rate=30)],
                      analyzers=[sti.PartnershipFormationAnalyzer()])
    else:
        sim = sti.Sim(n_agents=500, dur=5, rand_seed=1, diseases=[sti.HIV(init_prev=0.05)],
                      networks=[sti.MFNetwork()],
                      analyzers=[sti.PartnershipFormationAnalyzer()])
    sim.run(verbose=0)
    ana = sim.analyzers['partnershipformationanalyzer']
    rels = ana.results['relationships']['mfnetwork']
    assert len(rels) > 0, 'No edges recorded.'
    final_ti = ana._final_ti

    # See _MF_EDGES, above, to understand indicies into rel records
    finite = [rel for rel in rels if rel[3] is not None]
    assert len(finite) > 0, 'No expirations recorded — is super()._on_edge_dissolution missing?'
    assert all(rel[3] > rel[2] for rel in rels if rel[3] is not None), 'Found ti_expired <= ti_formed, minimum diff == 1.'
    assert all(rel[3] > rel[2] for rel in finite), 'Found ti_expired <= ti_formed.'
    assert all(rel[3] <= final_ti + 1 for rel in finite), 'Found ti_expired > final_ti + 1.'

    # Expirations recorded *during* the run (ti_expired <= final_ti), not all
    # deferred to the final-step sentinel. Catches normal breakups + death breakups.
    n_during_run = sum([1 for rel in finite if rel[3] <= final_ti])
    assert n_during_run > 0, 'No expirations recorded during the run (override or rel durs may be broken).'

    # Now check to ensure that the final state of the network is equivalently represented as the non-expired
    # network relationships, meaning, have a ti_expired == None .
    nw = sim.networks['mfnetwork']
    final_keys = set(zip(np.asarray(nw.edges.p1).tolist(),
                         np.asarray(nw.edges.p2).tolist(),
                         np.asarray(nw.edges.ti_formed).tolist()))
    none_keys = set((rel[0], rel[1], rel[2]) for rel in rels if rel[3] is None)
    assert len(none_keys) > 0, "There are no relationships left in network."
    assert none_keys == final_keys, (len(none_keys), len(final_keys))

    print(f'test_records_expired_end_to_end: {len(finite)}/{len(rels)} edges expired; '
          f'{len(none_keys)} active-at-end match the network. PASSED.')
    return sim


@sc.timer()
def test_deaths_intervals_consistent():
    """With deaths enabled, every removed edge (including death-driven, via the
    remove_uids override) is captured: no edge is left with ti_expired=None
    unless it is genuinely still active in the network at run end. Also checks
    intervals are valid and that expirations are actually recorded during the
    run (not all swept to the final sentinel) -- this is what regresses if the
    remove_uids override breaks. Guards against a removal path bypassing
    expiration recording."""
    sc.heading("Deaths: no edge stuck active; None-edges match the network active set")
    return test_records_expired_end_to_end(deaths=True)


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    # Injection-based getter tests (exact-value verification, no sim)
    test_get_unique_partners_active_per_agent()
    test_get_n_partnerships_formed_per_agent()
    test_get_n_partnerships_formed_binned()
    test_active_vs_formed()
    test_half_open_window_boundary()
    test_same_sex_network()
    test_multiple_networks()

    # Construction- and sim-level behavior tests
    test_duplicate_age_bins_rejected()
    test_network_missing_raises()
    test_network_nonsexual_warns_and_skips()
    test_network_unsupported_expiry_warns_and_skips()
    test_records_expired_end_to_end()
    test_deaths_intervals_consistent()

    # End-to-end integration (real multi-seed sim + plots/tables)
    test_partnership_formation_analyzer()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
