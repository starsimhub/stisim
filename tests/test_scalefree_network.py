"""Tests for ``sti.MSMScaleFreeNetwork`` and its FastSampler helper."""
import numpy as np
import pytest
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.networks._fast_sampler import FastSampler


# %% FastSampler unit tests

def test_fast_sampler_uniform_recovery():
    """Uniform weights → uniform samples (within tolerance)."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.ones(4))
    n_draws = 10_000
    counts = np.bincount([sampler.sample_index(rng) for _ in range(n_draws)], minlength=4)
    rel_freq = counts / n_draws
    assert np.all(np.abs(rel_freq - 0.25) < 0.02), f'non-uniform rel_freq: {rel_freq}'


def test_fast_sampler_weighted_recovery():
    """Weight ratio is preserved in the empirical sample."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.array([1.0, 3.0]))
    n_draws = 10_000
    counts = np.bincount([sampler.sample_index(rng) for _ in range(n_draws)], minlength=2)
    rel_freq = counts / n_draws
    # Expected: index 0 -> 0.25, index 1 -> 0.75
    assert abs(rel_freq[0] - 0.25) < 0.02
    assert abs(rel_freq[1] - 0.75) < 0.02


def test_fast_sampler_zero_mass_raises():
    """Empty or all-zero input raises on sample."""
    rng = np.random.default_rng(0)
    for w in [np.array([]), np.zeros(5)]:
        sampler = FastSampler(w)
        with pytest.raises(RuntimeError, match='no mass'):
            sampler.sample_index(rng)


def test_fast_sampler_handles_negative_weights():
    """Negative weights are clipped to zero (not raise)."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.array([-1.0, 2.0, -3.0, 4.0]))
    # Only indices 1 and 3 should ever be drawn
    drawn = {sampler.sample_index(rng) for _ in range(200)}
    assert drawn <= {1, 3}, f'drew clipped-out indices: {drawn}'


# %% MSMScaleFreeNetwork tests


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
    # ss.dur is stored as-is; converted to steps at init_pre.
    assert net._target_mean_dur_steps is None  # not yet converted


@sc.timer()
def test_init_pre_converts_durs_to_steps():
    """After sim init, dur params are converted to integer step counts."""
    net = sti.MSMScaleFreeNetwork()
    sim = _make_sim(net=net)
    sim.init()
    # sti.Sim deep-copies modules by default; pull the live network from the sim.
    sim_net = sim.networks[0]
    # default ss.years(2) and a default dt of 1/12 yr -> 24 steps
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
    pool = sim_net._get_pool()
    pool_uids = pool.uids
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
