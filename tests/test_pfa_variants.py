"""Smoke tests for pair-formation algorithm variants."""
import numpy as np
import starsim as ss
import stisim as sti
from stisim.pfa_variants import (
    MFNetwork_LSA,
    MFNetwork_SortPair,
    MFNetwork_DesiredAgeBucket,
    MFNetwork_SortBisect,
    MFNetwork_GreedyOldEnough,
)


def _match_once(net, n_agents=1_000, seed=0):
    """Init a sim and call match_pairs() directly. Returns (sim, p1, p2)."""
    sim = ss.Sim(n_agents=n_agents, networks=net, diseases='sis',
                 start='2000-01-01', stop='2001-01-01', rand_seed=seed)
    sim.init()
    net = sim.networks[0]  # sim deep-copies networks
    try:
        p1, p2 = net.match_pairs()
    except sti.networks.NoPartnersFound:
        empty = ss.uids(np.array([], dtype=np.int64))
        p1, p2 = empty, empty
    return sim, p1, p2


def _assert_valid_pairs(p1, p2, sim):
    assert len(p1) == len(p2)
    if len(p1) == 0:
        return
    assert (sim.people.male[p1]).all(), "p1 must be all male"
    assert (sim.people.female[p2]).all(), "p2 must be all female"


def test_lsa_variant_runs():
    net = MFNetwork_LSA()
    sim, p1, p2 = _match_once(net, n_agents=500)
    assert len(p1) > 0, "LSA variant should produce some edges"
    _assert_valid_pairs(p1, p2, sim)
    # LSA produces a strict 1-1 assignment: no duplicates.
    assert len(np.unique(p1)) == len(p1)
    assert len(np.unique(p2)) == len(p2)


def test_sortpair_variant_runs():
    net = MFNetwork_SortPair()
    sim, p1, p2 = _match_once(net, n_agents=500)
    assert len(p1) > 0
    _assert_valid_pairs(p1, p2, sim)
    # SortPair has no replacement -- p1 and p2 should be unique.
    assert len(np.unique(p1)) == len(p1)
    assert len(np.unique(p2)) == len(p2)


def test_desired_age_bucket_variant_runs():
    net = MFNetwork_DesiredAgeBucket()
    sim, p1, p2 = _match_once(net, n_agents=500)
    # Variant may post-filter aggressively; allow zero matches but if any, validate them.
    _assert_valid_pairs(p1, p2, sim)
    # No assertion on uniqueness: with sample-with-replacement, p1 can have duplicates
    # (capped by male concurrency).


def test_desired_age_bucket_post_filter_wired():
    """With all acceptance probabilities zero, the Bernoulli post-filter must drain output."""
    net = MFNetwork_DesiredAgeBucket(
        p_matched_stable=[0, 0, 0],
        p_mismatched_casual=[0, 0, 0],
    )
    sim, p1, p2 = _match_once(net, n_agents=500)
    assert len(p1) == 0, "all-zero acceptance should yield no matches"


def test_sortbisect_variant_runs():
    net = MFNetwork_SortBisect()
    sim, p1, p2 = _match_once(net, n_agents=500)
    assert len(p1) > 0
    _assert_valid_pairs(p1, p2, sim)
    # SortBisect (production) has no replacement.
    assert len(np.unique(p1)) == len(p1)
    assert len(np.unique(p2)) == len(p2)


def test_greedy_old_enough_runs():
    net = MFNetwork_GreedyOldEnough()
    sim, p1, p2 = _match_once(net, n_agents=500)
    _assert_valid_pairs(p1, p2, sim)
    if len(p1):
        # No replacement: each man appears at most once in p1.
        assert len(np.unique(p1)) == len(p1)
        assert len(np.unique(p2)) == len(p2)
