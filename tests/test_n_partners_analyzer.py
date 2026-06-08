"""
Unit tests for NPartnersAnalyzer (defined in
tests/devtests/devtest_network_diagnostics.py).

These tests inject controlled per-timestep edge records directly into an
analyzer instance (bypassing a full sim run) so that uniqueness, non-uniqueness,
duplicate handling, and the trailing-window boundary can be verified
deterministically for both male- and female-subject modes.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import sciris as sc

from stisim.analyzers import NPartnersAnalyzer

# NPartnersAnalyzer lives in the devtests diagnostics module; make it importable.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'devtests'))


verbose = False
do_plot = False
sc.options(interactive=False)


def _make_npartners_analyzer(window_months, records, male_uids, female_uids, final_ti):
    """Build an NPartnersAnalyzer with manually-injected state for unit testing.

    Each record is a 6-tuple of parallel sequences (one entry per edge active
    at that ``ti``)::

        (ti, male_uids, female_uids, formation_ti, age_male, age_female)

    where ``age_male`` / ``age_female`` are the partners' ages at formation (as
    recorded on the edge). This bypasses running a full sim so we can control
    duplicates, ages, and the window precisely.
    """
    ana = NPartnersAnalyzer(window_months=window_months)
    ana._records = [
        (ti,
         np.asarray(m, dtype=np.int64),
         np.asarray(f, dtype=np.int64),
         np.asarray(fti, dtype=np.int64),
         np.asarray(age_m, dtype=float),
         np.asarray(age_f, dtype=float))
        for (ti, m, f, fti, age_m, age_f) in records
    ]
    ana._final_male_uids = np.asarray(male_uids, dtype=np.int64)
    ana._final_female_uids = np.asarray(female_uids, dtype=np.int64)
    ana._final_ti = final_ti
    return ana


# Shared scenario for all three tests.
# window_months=3, final_ti=5 -> threshold = 5 - 3 + 1 = 3 (count ti in {3,4,5}).
# Males: 10, 11. Females: 20, 21, 22.
# Ages are recorded at formation; they vary by ti for the persisting edges but
# only the formation-step ages (where formation_ti == ti) are ever binned.
_SCENARIO = dict(
    window_months=3,
    male_uids=[10, 11],
    female_uids=[20, 21, 22],
    final_ti=5,
    records=[
        # (ti, males, females, formation_ti, age_male, age_female)
        # ti=2 is OUTSIDE the window and must be ignored entirely.
        (2, [10],             [20],             [2],          [30],             [22]),
        # ti=3: male 10 paired with female 20 TWICE (duplicate pair, same step)
        #       plus female 21; male 11 paired with female 22 (formed earlier, fti=1).
        (3, [10, 10, 10, 11], [20, 20, 21, 22], [3, 3, 3, 1], [31, 31, 31, 41], [22, 22, 18, 27]),
        # ti=4: the 10-20 and 11-22 edges persist (not newly formed).
        (4, [10, 11],         [20, 22],         [3, 1],       [32, 42],         [23, 28]),
        # ti=5: male 10 forms a new edge with female 22.
        (5, [10],             [22],             [5],          [33],             [29]),
    ],
)


@sc.timer()
def test_npartners_unique_partner_counts():
    """get_unique_partner_counts: unique opposite-sex partners per agent over the
    trailing window, for both male- and female-subject modes, with duplicate
    (male, female) appearances collapsed to a single unique partner."""
    sc.heading("Ensuring get_unique_partner_counts returns unique partner counts "
               "per male and per female, collapsing duplicate pairings")
    ana = _make_npartners_analyzer(**_SCENARIO)

    # Per male: unique females seen in window. Male 10 -> {20,21,22}=3 (20 seen
    # 3x but unique), male 11 -> {22}=1.
    male_counts = ana.get_unique_partner_counts(female=False)
    assert male_counts.tolist() == [3, 1], male_counts.tolist()

    # Per female: unique males seen in window. Female 20 -> {10}=1 (10 seen 3x
    # but unique), female 21 -> {10}=1, female 22 -> {11,10}=2.
    female_counts = ana.get_unique_partner_counts(female=True)
    assert female_counts.tolist() == [1, 1, 2], female_counts.tolist()

    print('test_npartners_unique_partner_counts: '
          f'male={male_counts.tolist()}, female={female_counts.tolist()}. PASSED.')
    return


@sc.timer()
def test_npartners_n_partnerships_formed():
    """get_n_partnerships_formed: total (non-unique) partnerships formed per
    agent within the trailing window, for both subject modes. Repeated
    partnerships with the same partner are each counted; only edges whose
    formation_ti falls in the window are counted (once, at formation)."""
    sc.heading("Ensuring get_n_partnerships_formed returns non-unique formation "
               "counts per male and per female within the trailing window")
    ana = _make_npartners_analyzer(**_SCENARIO)

    # Per male: formations with formation_ti in window. Male 10: at ti=3 three
    # new edges (10-20, 10-20, 10-21) + at ti=5 one (10-22) = 4 (the duplicate
    # 10-20 pair counts twice). Male 11: its only edge formed at ti=1 (outside
    # window) -> 0.
    male_formed = ana.get_n_partnerships_formed(female=False)
    assert male_formed.tolist() == [4, 0], male_formed.tolist()

    # Per female: female 20 -> 2 (both ti=3 formations), female 21 -> 1 (ti=3),
    # female 22 -> 1 (ti=5; its ti=1 formation is outside the window).
    female_formed = ana.get_n_partnerships_formed(female=True)
    assert female_formed.tolist() == [2, 1, 1], female_formed.tolist()

    print('test_npartners_n_partnerships_formed: '
          f'male={male_formed.tolist()}, female={female_formed.tolist()}. PASSED.')
    return


@sc.timer()
def test_npartners_n_partnerships_formed_binned():
    """get_n_partnerships_formed(age_bins=...): non-unique formations attributed
    to the subject's age *at formation*, with overlapping bins each counting a
    formation. Uses the same scenario; in-window formation ages are male
    [31, 31, 31, 33] and female [22, 22, 18, 29]."""
    sc.heading("Ensuring get_n_partnerships_formed bins formations by age at "
               "formation, with overlapping bins each counting")
    ana = _make_npartners_analyzer(**_SCENARIO)

    # Disjoint bins: 15-30 and 30-45.
    disjoint = ['15-30', '30-45']
    male_binned = ana.get_n_partnerships_formed(female=False, age_bins=disjoint)
    female_binned = ana.get_n_partnerships_formed(female=True, age_bins=disjoint)
    # Males formed at ages [31,31,31,33] -> all in 30-45: [0, 4].
    assert male_binned.tolist() == [0, 4], male_binned.tolist()
    # Females formed at ages [22,22,18,29] -> all in 15-30: [4, 0].
    assert female_binned.tolist() == [4, 0], female_binned.tolist()

    # Per-bin totals must match the (binless) overall formation totals when the
    # bins partition the observed ages without overlap.
    assert male_binned.sum() == ana.get_n_partnerships_formed(female=False).sum()
    assert female_binned.sum() == ana.get_n_partnerships_formed(female=True).sum()

    # Overlapping bins: 15-45 fully contains 30-45, so a 30-45 formation counts
    # in BOTH bins.
    overlapping = ['15-45', '30-45']
    male_overlap = ana.get_n_partnerships_formed(female=False, age_bins=overlapping)
    female_overlap = ana.get_n_partnerships_formed(female=True, age_bins=overlapping)
    # Males [31,31,31,33]: 15-45 -> 4, 30-45 -> 4.
    assert male_overlap.tolist() == [4, 4], male_overlap.tolist()
    # Females [22,22,18,29]: 15-45 -> 4, 30-45 -> 0 (all under 30).
    assert female_overlap.tolist() == [4, 0], female_overlap.tolist()

    print('test_npartners_n_partnerships_formed_binned: '
          f'disjoint male={male_binned.tolist()}, female={female_binned.tolist()}; '
          f'overlap male={male_overlap.tolist()}, female={female_overlap.tolist()}. PASSED.')
    return


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)
    timer = sc.timer()

    test_npartners_unique_partner_counts()
    test_npartners_n_partnerships_formed()
    test_npartners_n_partnerships_formed_binned()

    sc.heading("Total:")
    timer.toc()

    if do_plot:
        plt.show()
