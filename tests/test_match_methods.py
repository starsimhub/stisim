"""Tests for the match_method API on MFNetwork."""
import hashlib
import numpy as np
import pytest
import starsim as ss
import stisim as sti
from stisim.networks import matchers


def _make_sim(match_method=None, n_agents=2_000, seed=0):
    kw = {} if match_method is None else {'match_method': match_method}
    net = sti.MFNetwork(**kw)
    sim = ss.Sim(n_agents=n_agents, networks=net, diseases=sti.HIV(),
                 start='2000-01-01', stop='2001-01-01', rand_seed=seed)
    sim.init()
    return sim


def _fingerprint(p1, p2):
    """Stable hash over matched pairs for cross-version comparison."""
    a = np.asarray(p1, dtype=np.int64)
    b = np.asarray(p2, dtype=np.int64)
    raw = a.tobytes() + b'|' + b.tobytes()
    return hashlib.sha256(raw).hexdigest()


def test_sort_bisect_dispatch_is_reproducible():
    """sort_bisect dispatch is deterministic: same seed → same fingerprint."""
    sim_a = _make_sim(match_method='sort_bisect', seed=42)
    net_a = sim_a.networks[0]
    np.random.seed(0)
    p1_a, p2_a = net_a.match_pairs()

    sim_b = _make_sim(match_method='sort_bisect', seed=42)
    net_b = sim_b.networks[0]
    np.random.seed(0)
    p1_b, p2_b = net_b.match_pairs()

    assert _fingerprint(p1_a, p2_a) == _fingerprint(p1_b, p2_b)


@pytest.mark.parametrize('method', sorted(matchers.MATCHERS.keys()))
def test_each_matcher_produces_valid_pairs(method):
    """Smoke test: each registered matcher runs and produces valid (p1, p2)."""
    seed = abs(hash(method)) % 1000
    sim = _make_sim(match_method=method, n_agents=500, seed=seed)
    net = sim.networks[0]
    try:
        p1, p2 = net.match_pairs()
    except sti.networks.NoPartnersFound:
        # desired_age_bucket can post-filter to empty at small n; that's fine.
        return
    assert len(p1) == len(p2)
    if len(p1):
        assert sim.people.male[p1].all(), 'p1 must be all male'
        assert sim.people.female[p2].all(), 'p2 must be all female'


def test_callable_match_method_is_invoked():
    """Passing a callable as match_method bypasses the registry."""
    sentinel = {'called': False}

    def custom(net):
        sentinel['called'] = True
        empty = ss.uids(np.array([], dtype=np.int64))
        return empty, empty

    sim = _make_sim(match_method=custom, n_agents=200, seed=0)
    net = sim.networks[0]
    try:
        net.match_pairs()
    except sti.networks.NoPartnersFound:
        pass
    assert sentinel['called'], 'custom matcher must be called'


def test_unknown_string_raises_keyerror():
    """Unknown method strings raise KeyError."""
    sim = _make_sim(match_method='no_such_method', n_agents=200, seed=0)
    net = sim.networks[0]
    with pytest.raises(KeyError):
        net.match_pairs()


def test_default_match_method_is_kdtree_nn():
    """A bare MFNetwork() uses kdtree_nn by default."""
    net = sti.MFNetwork()
    assert net.pars.match_method == 'kdtree_nn'
