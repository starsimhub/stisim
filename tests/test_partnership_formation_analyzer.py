"""
Tests for stisim.PartnershipFormationAnalyzer.

Runs a small STI simulation, prints tabular partnership-formation counts per
gender / age bin / timestep, and produces two plots:
  1. Total male vs female new partnerships per timestep.
  2. Stacked per-age-bin subplots (youngest at bottom, oldest at top), each
     with male/female lines.
"""

import numpy as np
import pandas as pd
import pylab as pl
import starsim as ss
import stisim as sti


# Default test configuration (configurable via function arguments)
DEFAULT_N_AGENTS  = 1000
DEFAULT_DUR_YEARS = 5
DEFAULT_AGE_BINS  = [f'{lo}-{lo+5}' for lo in range(0, 75, 5)]  # 0-5 ... 70-75

# NOTE: the instructions call for dt=1/12 (monthly), but sti.Sim's route_pars
# broadcasts ``dt`` to ['sim', 'sti', 'connector'] which then trips either the
# disease-instance XOR check or the unused-connector-par check. Until sti.Sim
# handles universal pars cleanly, we let the sim use its default dt.

# Plot colors per the instructions
MALE_COLOR   = 'blue'
FEMALE_COLOR = 'red'


ANALYZER_NAME = 'partnershipformationanalyzer'


def make_sim(n_agents=DEFAULT_N_AGENTS, dur=DEFAULT_DUR_YEARS, age_bins=None):
    """Construct a sim with a StructuredSexual network and the analyzer.

    Returns the sim only — fetch the analyzer post-init via
    ``sim.analyzers[ANALYZER_NAME]`` because sti.Sim deep-copies the analyzer
    during init, leaving any pre-init reference stale.
    """
    age_bins = list(age_bins) if age_bins is not None else list(DEFAULT_AGE_BINS)
    sim = sti.Sim(
        n_agents=n_agents,
        dur=dur,
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
    """Print one wide table per sex: rows = age bins, columns = timesteps."""
    for sex in ('f', 'm'):
        label = 'FEMALE' if sex == 'f' else 'MALE'
        sub = df[df.sex == sex]
        wide = sub.pivot(index='age_bin', columns='ti', values='count')
        wide = wide.reindex(analyzer.age_bins)  # keep age order
        print(f'\n=== New partnerships per timestep ({label}) ===')
        with pd.option_context('display.max_columns', None,
                               'display.width', 200,
                               'display.max_rows', None):
            print(wide)


def plot_totals(df, sim):
    """Total new partnerships per timestep, one line per sex."""
    yearvec = sim.t.yearvec
    n_ti = len(yearvec)
    male   = (df[df.sex == 'm'].groupby('ti')['count'].sum()
                .reindex(range(n_ti)).fillna(0).values)
    female = (df[df.sex == 'f'].groupby('ti')['count'].sum()
                .reindex(range(n_ti)).fillna(0).values)

    fig, ax = pl.subplots(figsize=(9, 5))
    ax.plot(yearvec, male,   color=MALE_COLOR,   label='Male')
    ax.plot(yearvec, female, color=FEMALE_COLOR, label='Female')
    ax.set_xlabel('Year')
    ax.set_ylabel('Number of partnerships formed')
    ax.set_title('Total partnerships formed per timestep')
    ax.legend()
    pl.tight_layout()
    return fig


def plot_by_age_bin(df, analyzer, sim):
    """Stacked subplots (one per age bin), youngest at bottom, oldest at top."""
    yearvec = sim.t.yearvec
    n_bins = len(analyzer.age_bins)
    fig, axes = pl.subplots(n_bins, 1, figsize=(9, max(1.2 * n_bins, 4)),
                            sharex=True)
    if n_bins == 1:
        axes = [axes]

    for i, bin_str in enumerate(analyzer.age_bins):
        # axes[0] is the TOP row. Youngest -> bottom -> last axis.
        ax = axes[n_bins - 1 - i]
        sub = df[df.age_bin == bin_str].sort_values('ti')
        male   = sub[sub.sex == 'm']['count'].values
        female = sub[sub.sex == 'f']['count'].values
        ax.plot(yearvec, male,   color=MALE_COLOR,   label='Male')
        ax.plot(yearvec, female, color=FEMALE_COLOR, label='Female')
        ax.set_ylabel(bin_str, fontsize=9, rotation=0, ha='right', va='center')

    axes[-1].set_xlabel('Year')
    axes[0].set_title('Partnerships formed per timestep by age bin '
                      '(youngest at bottom)')
    axes[0].legend(loc='upper right', fontsize=8)
    pl.tight_layout()
    return fig


def test_partnership_formation_analyzer(n_agents=DEFAULT_N_AGENTS,
                                        dur=DEFAULT_DUR_YEARS,
                                        age_bins=None,
                                        show_plots=False):
    """End-to-end test: run sim, check counts, print tables, build plots."""
    sim = make_sim(n_agents=n_agents, dur=dur, age_bins=age_bins)
    sim.run()
    # sti.Sim deep-copies analyzer instances during init, so fetch the live
    # post-run analyzer from sim.analyzers (lookup is by ss.Module.name).
    analyzer = sim.analyzers['partnershipformationanalyzer']

    df = build_dataframe(sim, analyzer)

    # ----- Sanity checks -----
    total = int(df['count'].sum())
    assert total > 0, 'Expected at least some new partnerships in a 5-year sim.'

    # StructuredSexual is male-female only. When both partners fall inside the
    # age bins, total male and total female counts must match exactly. Allow a
    # small tolerance for partnerships where one partner ages above the top bin.
    male_total   = int(df[df.sex == 'm']['count'].sum())
    female_total = int(df[df.sex == 'f']['count'].sum())
    max_total = max(male_total, female_total, 1)
    diff_share = abs(male_total - female_total) / max_total
    assert diff_share < 0.10, (
        f'Male/female partnership totals diverge by more than 10% '
        f'(m={male_total}, f={female_total}). For a male-female-only network '
        f'they should be approximately equal.'
    )

    # No counts should appear in pre-debut age bins (0-5, 5-10) — debut age is
    # well above 10 for the StructuredSexual defaults.
    pre_debut_count = int(df[df.age_bin.isin(['0-5', '5-10'])]['count'].sum())
    assert pre_debut_count == 0, (
        f'Found {pre_debut_count} partnerships in pre-debut age bins (0-5, 5-10).'
    )

    # First timestep typically has no new partnerships (sim just initialized);
    # but later timesteps should accumulate counts.
    per_ti_total = df.groupby('ti')['count'].sum().values
    assert per_ti_total[1:].sum() > 0, 'No new partnerships after the first timestep.'

    # ----- Reporting -----
    print(f'\nSim: n_agents={n_agents}, dur={dur} years, '
          f'{len(sim.t.yearvec)} timesteps (default sim dt)')
    print(f'Total new partnerships recorded (male+female partner counts): {total}')
    print(f'  Male contributions:   {male_total}')
    print(f'  Female contributions: {female_total}')

    print_tables(df, analyzer)
    fig_totals = plot_totals(df, sim)
    fig_by_bin = plot_by_age_bin(df, analyzer, sim)

    if show_plots:
        pl.show()
    else:
        pl.close(fig_totals)
        pl.close(fig_by_bin)

    return sim, analyzer, df


def _first_bin_for_age(ana, age):
    """Return the index of the first age bin containing ``age``, or -1 if none."""
    for i, (lo, hi) in enumerate(ana.age_ranges):
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
    ana = sim.analyzers['partnershipformationanalyzer']
    ppl = sim.people

    # Find one female and one male UID whose age falls inside some bin
    female_uid = None
    male_uid = None
    for uid in range(n_agents):
        if _first_bin_for_age(ana, float(ppl.age[uid])) < 0:
            continue
        if bool(ppl.female[uid]) and female_uid is None:
            female_uid = uid
        elif (not bool(ppl.female[uid])) and male_uid is None:
            male_uid = uid
        if female_uid is not None and male_uid is not None:
            break
    assert female_uid is not None and male_uid is not None, (
        'Could not find both a female and male UID in-range; raise n_agents.'
    )

    # female_uid appears 3 times, male_uid appears 1 time
    uid_array = np.array([female_uid, male_uid, female_uid, female_uid])
    f_bin = _first_bin_for_age(ana, float(ppl.age[female_uid]))
    m_bin = _first_bin_for_age(ana, float(ppl.age[male_uid]))

    f_counts = ana._count_by_age_bin(uid_array, female=True)
    m_counts = ana._count_by_age_bin(uid_array, female=False)

    assert f_counts[f_bin] == 3, (
        f'Female UID appeared 3 times but bin {f_bin} count is {f_counts[f_bin]} '
        f'(full counts: {f_counts.tolist()}).'
    )
    assert m_counts[m_bin] == 1, (
        f'Male UID appeared 1 time but bin {m_bin} count is {m_counts[m_bin]} '
        f'(full counts: {m_counts.tolist()}).'
    )
    # No spurious counts elsewhere (default bins are non-overlapping)
    assert int(f_counts.sum()) == 3, f'Stray female counts: {f_counts.tolist()}'
    assert int(m_counts.sum()) == 1, f'Stray male counts: {m_counts.tolist()}'

    # Empty input → all-zero counts
    empty = np.empty(0, dtype=int)
    assert int(ana._count_by_age_bin(empty, female=True).sum()) == 0
    assert int(ana._count_by_age_bin(empty, female=False).sum()) == 0

    print(f'test_count_by_age_bin_duplicates: female bin {f_bin} = 3, '
          f'male bin {m_bin} = 1, empty input = 0. PASSED.')
    return


def test_count_by_age_bin_overlapping(n_agents=200):
    """Regression: with overlapping age bins, each agent is counted in every
    bin whose range contains their age (not just one).
    """
    overlapping_bins = ['15-50', '40-50']  # 40-49 falls in BOTH
    sim = sti.Sim(
        n_agents=n_agents, dur=1,
        diseases=[sti.HIV(init_prev=0.05)],
        networks=[sti.StructuredSexual()],
        analyzers=[sti.PartnershipFormationAnalyzer(age_bins=overlapping_bins)],
    )
    sim.init()
    ana = sim.analyzers['partnershipformationanalyzer']
    ppl = sim.people
    assert ana.age_bins == overlapping_bins
    assert ana.age_ranges == [(15.0, 50.0), (40.0, 50.0)]

    # Find a female agent aged 40-49 (in BOTH bins)
    both_bins_female = None
    for uid in range(n_agents):
        age = float(ppl.age[uid])
        if bool(ppl.female[uid]) and 40 <= age < 50:
            both_bins_female = uid
            break
    assert both_bins_female is not None, (
        'Could not find a female aged 40-49; raise n_agents.'
    )

    # Single occurrence: expect count of 1 in BOTH bins
    counts = ana._count_by_age_bin(np.array([both_bins_female]), female=True)
    assert counts.tolist() == [1, 1], (
        f'Overlapping-bins fix: agent aged 40-49 should count in both bins, '
        f'got {counts.tolist()} (expected [1, 1]).'
    )

    # Repeated occurrences should still accumulate independently in each bin
    counts_repeated = ana._count_by_age_bin(
        np.array([both_bins_female, both_bins_female, both_bins_female]), female=True)
    assert counts_repeated.tolist() == [3, 3], (
        f'Three occurrences should give [3, 3], got {counts_repeated.tolist()}.'
    )

    # A female aged 15-39 should only land in the wide bin
    only_wide_female = None
    for uid in range(n_agents):
        age = float(ppl.age[uid])
        if bool(ppl.female[uid]) and 15 <= age < 40:
            only_wide_female = uid
            break
    if only_wide_female is not None:
        counts_one = ana._count_by_age_bin(np.array([only_wide_female]), female=True)
        assert counts_one.tolist() == [1, 0], (
            f'Female aged 15-39 should land only in 15-50, got {counts_one.tolist()}.'
        )

    print(f'test_count_by_age_bin_overlapping: agent aged 40-49 counted in both bins, '
          f'repeated occurrences accumulate. PASSED.')
    return


if __name__ == '__main__':
    test_count_by_age_bin_duplicates()
    test_count_by_age_bin_overlapping()
    test_partnership_formation_analyzer(show_plots=True)