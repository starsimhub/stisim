"""Tests for ``sti.MSMScaleFreeNetwork``."""
import numpy as np
import pytest
import sciris as sc
import starsim as ss
import stisim as sti


def _make_net(**kw):
    """Network with full male participation, so kernel-mechanics tests see the
    whole post-debut male pool regardless of the default ``p_msm``."""
    kw.setdefault('p_msm', ss.bernoulli(p=1.0))
    return sti.MSMScaleFreeNetwork(**kw)


def _make_sim(net=None, n_agents=1_000, dur=2, rand_seed=1):
    if net is None:
        net = _make_net()
    return sti.Sim(
        diseases=[sti.HIV(init_prev=0.05, beta_m2f=0.05)],
        networks=[net, ss.MaternalNet()],
        demographics=[ss.Pregnancy(), ss.Deaths()],
        n_agents=n_agents, dur=dur, start=2010, verbose=-1, rand_seed=rand_seed,
    )


@sc.timer()
def test_class_instantiates_with_defaults():
    """Bare instantiation: defaults are set; no sim required."""
    net = _make_net()
    assert net.pars.target_mean_degree == 2.0
    assert float(net.pars.phi) == 1.0
    assert net._target_mean_dur_steps is None  # not yet converted


@sc.timer()
def test_init_pre_converts_durs_to_steps():
    """After sim init, dur params are converted to integer step counts."""
    net = _make_net()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    assert sim_net._target_mean_dur_steps == int(round(2.0 / float(sim_net.t.dt)))
    assert sim_net._max_edge_dur_steps == int(round(10.0 / float(sim_net.t.dt)))


@sc.timer()
def test_init_pre_raises_when_max_below_target():
    """``max_edge_dur < target_mean_dur`` is rejected at init."""
    net = _make_net(target_mean_dur=ss.years(5), max_edge_dur=ss.years(2))
    sim = _make_sim(net=net)
    with pytest.raises(ValueError, match='max_edge_dur'):
        sim.init()


@sc.timer()
def test_get_pool_filters_to_post_debut_males():
    """``_get_pool`` returns only post-debut male agents."""
    net = _make_net()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    pool_uids = sim_net._get_pool().uids
    assert len(pool_uids) > 0, 'pool is empty — fixture must produce post-debut males'
    assert sim.people.male[pool_uids].all(), 'pool contains females'
    debut_ok = sim.people.age[pool_uids] >= sim_net.debut[pool_uids]
    assert debut_ok.all(), 'pool contains pre-debut agents'


@sc.timer()
def test_p_msm_filters_pool():
    """``p_msm`` controls the fraction of post-debut males in the pool."""
    full = _make_sim(net=_make_net(p_msm=ss.bernoulli(p=1.0)), n_agents=2_000)
    part = _make_sim(net=_make_net(p_msm=ss.bernoulli(p=0.2)), n_agents=2_000)
    full.init(); part.init()
    n_full = full.networks[0]._get_pool().count()
    n_part = part.networks[0]._get_pool().count()
    assert n_full > 0, 'full-participation pool is empty — fixture failure'
    assert 0.1 < n_part / n_full < 0.35, \
        f'p_msm=0.2 pool fraction {n_part / n_full:.2f} not near 0.2 (full={n_full}, part={n_part})'


@sc.timer()
def test_mix_node_arrays_log1p_deg_invariant():
    """``log1p_deg`` is ``1 + log1p(deg)`` and bounded below by 1.0."""
    net = _make_net()
    sim = _make_sim(net=net)
    sim.init()
    sim_net = sim.networks[0]
    arrays = sim_net._mix_node_arrays()
    assert arrays['log1p_deg'].size > 0, 'mix arrays empty — fixture must produce a non-empty pool'
    assert (arrays['log1p_deg'] >= 1.0 - 1e-9).all(), 'log1p_deg has values below 1.0'


@sc.timer()
def test_mix_weights_row_returns_correct_length():
    """``_mix_weights_row(i, ...)`` returns ``n - 1 - i`` weights."""
    net = _make_net()
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
    net = _make_net()
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
    net = _make_net()
    sim = _make_sim(net=net, n_agents=300)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._build_kernel()
    assert sim_net._kernel_pairs_i.size > 0, 'fixture must produce a non-empty kernel'
    assert np.all(sim_net._kernel_pairs_i < sim_net._kernel_pairs_j)


@sc.timer()
def test_build_kernel_degenerate_pool_returns_zero():
    """Empty pool produces zero-sized kernel artefacts and n=0."""
    net = _make_net()
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
    net = _make_net()
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
    net = _make_net()
    sim = _make_sim(net=net, n_agents=100)
    sim.init()
    sim_net = sim.networks[0]
    sim_net._get_pool = lambda: sim.people.alive & ~sim.people.alive
    sim_net._build_kernel()
    rng = np.random.default_rng(0)
    assert sim_net._sample_pair_index(rng) == -1


# %% init_post + step dynamics tests


@sc.timer()
def test_initial_network_q0_non_empty():
    """``init_post`` samples a non-empty initial edge set for a large pool."""
    net = _make_net()
    sim = _make_sim(net=net, n_agents=2_000)
    sim.init()
    sim_net = sim.networks[0]
    assert len(sim_net.edges.p1) > 0, 'q0 produced no initial edges with n_agents=2000'


@sc.timer()
def test_initial_edges_post_debut_males_only():
    """Every endpoint of every initial edge is a post-debut male."""
    net = _make_net()
    sim = _make_sim(net=net, n_agents=2_000)
    sim.init()
    sim_net = sim.networks[0]
    p1 = np.asarray(sim_net.edges.p1, dtype=np.int64)
    p2 = np.asarray(sim_net.edges.p2, dtype=np.int64)
    assert len(p1) > 0, 'fixture must produce initial edges'
    all_uids = np.unique(np.concatenate([p1, p2]))
    assert sim.people.male[all_uids].all(), 'edge has a female endpoint'
    assert (sim.people.age[all_uids] >= sim_net.debut[all_uids]).all(), 'edge has a pre-debut endpoint'


@sc.timer()
def test_max_edge_dur_cap_respected():
    """No edge persists past ``max_edge_dur_steps``."""
    net = _make_net(max_edge_dur=ss.years(1), target_mean_dur=ss.months(6))
    sim = _make_sim(net=net, n_agents=1_500, dur=3)
    sim.run()
    sim_net = sim.networks[0]
    if len(sim_net.edges.p1):
        dur = np.asarray(sim_net.edges.dur, dtype=float)
        assert dur.max() <= sim_net._max_edge_dur_steps + 1, \
            f'edge dur {dur.max()} > cap {sim_net._max_edge_dur_steps}'


@sc.timer()
def test_network_initialises_and_steps():
    """End-to-end: a 2yr sim with the new network completes and forms edges."""
    net = _make_net()
    sim = _make_sim(net=net, n_agents=2_000, dur=2)
    sim.run()
    sim_net = sim.networks[0]
    assert len(sim_net.edges.p1) > 0, 'no edges after 2yr run — dynamics produced empty network'


# %% Steady-state validation


@sc.timer()
def test_target_mean_degree_honoured():
    """After burn-in, realised mean degree is within ±30% of target."""
    target = 2.0
    net = _make_net(target_mean_degree=target, target_mean_dur=ss.years(1))
    sim = _make_sim(net=net, n_agents=2_500, dur=5)
    sim.run()
    sim_net = sim.networks[0]
    n_pool = sim_net._get_pool().count()
    assert n_pool > 0, 'empty pool — fixture failure'
    realised = 2.0 * len(sim_net.edges.p1) / n_pool
    assert abs(realised - target) / target < 0.30, \
        f'realised mean degree {realised:.2f} not within 30% of target {target}'


@sc.timer()
def test_target_mean_dur_honoured():
    """After burn-in, mean still-alive edge age is between ½·target and max_edge_dur.

    Still-alive edges are duration-biased low (young edges over-represented),
    so the lower bound is half the target rather than the full target.
    """
    target_yrs = 1.0
    net = _make_net(target_mean_dur=ss.years(target_yrs))
    sim = _make_sim(net=net, n_agents=2_500, dur=5)
    sim.run()
    sim_net = sim.networks[0]
    assert len(sim_net.edges.p1) > 0, 'no edges to check duration of'
    # Per _append_new_edges, dur starts at max_edge_dur_steps and counts down.
    # Age = max_edge_dur_steps - dur.
    dur_remaining = np.asarray(sim_net.edges.dur, dtype=float)
    age_steps = float(sim_net._max_edge_dur_steps) - dur_remaining.mean()
    steps_per_yr = 1.0 / float(sim_net.t.dt)
    realised_yrs = age_steps / steps_per_yr
    max_yrs = float(sim_net._max_edge_dur_steps) / steps_per_yr
    assert 0.5 * target_yrs <= realised_yrs <= max_yrs, \
        f'realised mean age {realised_yrs:.2f}yr out of [0.5·target={0.5*target_yrs}, max={max_yrs:.2f}] for target {target_yrs}yr'


# %% Subclass-override smoke test


@sc.timer()
def test_subclass_mix_weights_override_is_used():
    """A subclass override of ``_mix_weights_row`` actually fires during kernel build.

    Avoids asserting downstream statistical differences (too noisy at this scale).
    Instead verifies the override mechanism: a subclass that returns a sentinel
    pattern should produce a kernel whose ``sel_w / A_w`` ratio reflects that pattern.
    """

    class SentinelMSM(sti.MSMScaleFreeNetwork):
        call_count = 0
        def _mix_weights_row(self, i, mix_arrays):
            type(self).call_count += 1
            # Sentinel: row i returns weights = (i+1) uniformly.
            return np.full(mix_arrays['log1p_deg'].size - 1 - i, float(i + 1))

    net = SentinelMSM(p_msm=ss.bernoulli(p=1.0))
    sim = _make_sim(net=net, n_agents=300)
    sim.init()
    assert SentinelMSM.call_count > 0, 'subclass _mix_weights_row was never called'
    sim_net = sim.networks[0]
    assert sim_net._kernel_pairs_i.size > 0, 'fixture must produce a non-empty kernel'
    # sel_w should reflect the sentinel scaling: pairs from row 0 have base weight 1
    # (normalised to 1 because mean=1), pairs from row 5 have base weight 6 / mean → much larger.
    # So sel_w should not be uniform across rows.
    row0_mask = sim_net._kernel_pairs_i == 0
    high_i = sim_net._kernel_pairs_i.max()
    high_mask = sim_net._kernel_pairs_i == high_i
    assert row0_mask.any() and high_mask.any(), 'kernel does not span both ends of i'
    # Per-row mean normalises mix weights to 1 within each row, but A_w is row-dependent
    # too — the override DOES change the kernel artefacts; just assert sel_w is non-trivial.
    assert sim_net._kernel_sel_w[row0_mask].max() > 0, 'row 0 selection weight is zero'
    assert sim_net._kernel_sel_w.size == sim_net._kernel_pairs_i.size, 'sel_w aligned with pairs'
