"""Tests for ``sti.MSMScaleFreeNetwork``."""
import numpy as np
import pytest
import sciris as sc
import starsim as ss
import stisim as sti


def _make_sim(net=None, n_agents=1_000, dur=2, rand_seed=1):
    if net is None:
        net = sti.MSMScaleFreeNetwork()
    return sti.Sim(
        diseases=[sti.HIV(init_prev=0.05, beta_m2f=0.05)],
        networks=[net, ss.MaternalNet()],
        demographics=[ss.Pregnancy(), ss.Deaths()],
        n_agents=n_agents, dur=dur, start=2010, verbose=-1, rand_seed=rand_seed,
    )


@sc.timer()
def test_class_instantiates_with_defaults():
    """Bare instantiation: defaults are set; no sim required."""
    net = sti.MSMScaleFreeNetwork()
    assert net.pars.target_mean_degree == 2.0
    assert float(net.pars.phi) == 1.0
    assert net._target_mean_dur_steps is None  # not yet converted


@sc.timer()
def test_init_pre_converts_durs_to_steps():
    """After sim init, dur params are converted to integer step counts."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    assert sim_net._target_mean_dur_steps == int(round(2.0 / float(sim_net.t.dt)))
    assert sim_net._max_edge_dur_steps == int(round(10.0 / float(sim_net.t.dt)))


@sc.timer()
def test_init_pre_raises_when_max_below_target():
    """``max_edge_dur < target_mean_dur`` is rejected at init."""
    net = sti.MSMScaleFreeNetwork(target_mean_dur=ss.years(5), max_edge_dur=ss.years(2))
    sim = _make_sim(net=net)
    with pytest.raises(ValueError, match='max_edge_dur'):
        sim.init()


@sc.timer()
def test_get_pool_filters_to_post_debut_males():
    """``_get_pool`` returns only post-debut male agents."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    pool_uids = sim_net._get_pool().uids
    assert len(pool_uids) > 0, 'pool is empty — fixture must produce post-debut males'
    assert sim.people.male[pool_uids].all(), 'pool contains females'
    debut_ok = sim.people.age[pool_uids] >= sim_net.debut[pool_uids]
    assert debut_ok.all(), 'pool contains pre-debut agents'


@sc.timer()
def test_mix_node_arrays_zero_degree_at_init():
    """Before any edges form, every agent has log1p_deg == 1."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    arrays = sim_net._mix_node_arrays()
    assert arrays['log1p_deg'].size > 0, 'mix arrays empty — fixture must produce a non-empty pool'
    assert np.allclose(arrays['log1p_deg'], 1.0)


@sc.timer()
def test_mix_weights_row_returns_correct_length():
    """``_mix_weights_row(i, ...)`` returns ``n - 1 - i`` weights."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    arrays = sim_net._mix_node_arrays()
    n = arrays['log1p_deg'].size
    assert n >= 3, f'fixture pool size {n} too small for shape test (need >= 3)'
    for i in range(min(n - 1, 3)):
        w = sim_net._mix_weights_row(i, arrays)
        assert w.shape == (n - 1 - i,), f'i={i}: got {w.shape}, expected ({n - 1 - i},)'
        assert np.all(w >= 0.0), 'mix weights must be non-negative'


@sc.timer()
def test_build_kernel_populates_artefacts():
    """``_build_kernel`` populates all the caches after init."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net, n_agents=300)
    sim.init()
    sim_net = sim.networks[0]
    n = sim_net._build_kernel()
    assert n >= 2, f'fixture must produce a non-degenerate pool (got n={n})'
    assert sim_net._kernel_pairs_i.size == sim_net._kernel_pairs_j.size
    assert sim_net._kernel_pairs_i.size == sim_net._kernel_A_w.size
    assert sim_net._kernel_pairs_i.size == sim_net._kernel_sel_w.size
    assert sim_net._kernel_hat_Ra > 0
    assert 0 < sim_net._kernel_q0 < 1


@sc.timer()
def test_build_kernel_pair_indices_upper_triangle():
    """``pairs_i < pairs_j`` for every entry (upper triangle only)."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net, n_agents=300)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._build_kernel()
    assert sim_net._kernel_pairs_i.size > 0, 'fixture must produce a non-empty kernel'
    assert np.all(sim_net._kernel_pairs_i < sim_net._kernel_pairs_j)


@sc.timer()
def test_build_kernel_degenerate_pool_returns_zero():
    """Empty pool produces zero-sized kernel artefacts and n=0."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net, n_agents=100)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._get_pool = lambda: sim.people.alive & ~sim.people.alive
    n = sim_net._build_kernel()
    assert n == 0
    assert sim_net._kernel_pairs_i.size == 0


@sc.timer()
def test_sample_pair_index_in_range_and_weighted():
    """``_sample_pair_index`` returns valid indices weighted by sel_w."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net, n_agents=300)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._build_kernel()
    m = sim_net._kernel_sel_w.size
    assert m > 0, 'fixture must produce a non-empty selection cdf'
    rng = np.random.default_rng(0)
    draws = np.array([sim_net._sample_pair_index(rng) for _ in range(500)])
    assert (draws >= 0).all() and (draws < m).all()


@sc.timer()
def test_sample_pair_index_returns_neg_on_empty_kernel():
    """Empty kernel → sentinel -1 from ``_sample_pair_index``."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net, n_agents=100)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._get_pool = lambda: sim.people.alive & ~sim.people.alive
    sim_net._build_kernel()
    rng = np.random.default_rng(0)
    assert sim_net._sample_pair_index(rng) == -1
