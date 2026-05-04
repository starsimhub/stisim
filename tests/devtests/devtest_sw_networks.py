"""
Integration tests for the MF + SW network split (issue #307).

Demonstrates the network design introduced in v1.5.4:
- ``MFNetwork`` — heterosexual, risk-group structured
- ``SWNetwork`` — sex work; FSW–client edges with per-agent age window
- ``StructuredSexual`` — the standard MF + SW combined network most users
  will reach for (single edge list with all four edge types)

Records the kind of network diagnostics shown in the eswatini
``plot_fig3_network.py`` figure (mean lifetime partners by sex, partnership
durations by type, risk-group composition, SW shares) to verify network
dynamics behave as expected.

Also covers all ``sti.Sim`` and ``hivsim.Sim`` construction paths, plus
both ``hivsim.demo`` examples.

Run from the repo root:
    pytest tests/devtests/devtest_sw_networks.py -v
or as a script:
    python tests/devtests/devtest_sw_networks.py
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
import hivsim


# Common knobs — small/short for fast iteration. ``copy_inputs=True`` (the
# sti.Sim default) means we can share these list literals across sims.
N_AGENTS = 1000
DUR = 10
SEED = 1
DEMOG = [ss.Pregnancy(), ss.Deaths()]


# ---------- Diagnostic helpers ----------

def network_stats(nw, people):
    """Compute network diagnostics on a single network module.

    Returns a dict with: mean lifetime partners (by sex), mean intended
    relationship duration (by edge type, in timesteps), risk-group shares
    (if present), and SW state counts (if present).
    """
    out = sc.objdict()
    alive = people.alive.values

    if hasattr(nw, 'lifetime_partners'):
        lp = nw.lifetime_partners.values
        out.lifetime_partners_f_mean = float(lp[alive & people.female.values].mean())
        out.lifetime_partners_m_mean = float(lp[alive & people.male.values].mean())

    for key in ('stable', 'casual', 'onetime', 'sw'):
        attr = f'lifetime_{key}_partners'
        if hasattr(nw, attr):
            v = getattr(nw, attr).values[alive]
            out[f'{attr}_mean'] = float(v.mean())

    if hasattr(nw, 'risk_group'):
        rg = nw.risk_group.values[alive]
        out.risk_group_share = {f'rg{i}': float((rg == i).mean()) for i in range(3)}

    if getattr(nw, 'relationship_durs', None):
        durs_by_type = sc.ddict(list)
        type_lookup = {v: k for k, v in nw.edge_types.items()}
        for entries in nw.relationship_durs.values():
            for e in entries:
                durs_by_type[type_lookup.get(e['edge_type'], '?')].append(e['dur'])
        for k, vals in durs_by_type.items():
            out[f'mean_dur_{k}'] = float(np.mean(vals)) if vals else 0.0

    if hasattr(nw, 'ever_fsw'):
        out.ever_fsw_count = int(nw.ever_fsw.sum())
        out.ever_client_count = int(nw.ever_client.sum())
        out.current_fsw_count = int(nw.fsw.sum())
        out.current_client_count = int(nw.client.sum())

    if hasattr(nw.edges, 'edge_type') and nw.edge_types:
        for k, v in nw.edge_types.items():
            out[f'edges_{k}'] = int((nw.edges.edge_type == v).sum())

    return out


def print_summary(label, stats):
    print(f'\n--- {label} ---')
    for k, v in stats.items():
        if isinstance(v, dict):
            print(f'  {k}: {v}')
        elif isinstance(v, float):
            print(f'  {k}: {v:.3f}')
        else:
            print(f'  {k}: {v}')


# ---------- 1. MFNetwork only via sti.Sim ----------

def test_mf_only():
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sti.MFNetwork()],
        demographics=DEMOG,
        dur=DUR,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    nw = sim.networks.mfnetwork
    s = network_stats(nw, sim.people)
    print_summary('MFNetwork only', s)

    # Observed values @ N_AGENTS=1000, DUR=10, SEED=1: f_mean≈1.58, m_mean≈1.55,
    # stable≈0.55, casual≈0.37, onetime≈0.65, mean_dur_stable≈1020 months,
    # mean_dur_casual≈15 months, edges_stable≈270.
    assert 1.0 < s.lifetime_partners_f_mean < 2.5
    assert 1.0 < s.lifetime_partners_m_mean < 2.5
    assert 0.3 < s.lifetime_stable_partners_mean < 1.0
    assert 0.2 < s.lifetime_casual_partners_mean < 0.7
    assert 0.3 < s.lifetime_onetime_partners_mean < 1.2
    # Stable partnerships are long-lived (risk-group-0 cap = 100y); casual ~years
    assert 600 < s.mean_dur_stable < 1400
    assert 8 < s.mean_dur_casual < 30
    # Risk-group shares sum to 1; default is f0=0.85, f2=0.01
    assert abs(sum(s.risk_group_share.values()) - 1.0) < 1e-6
    assert s.risk_group_share['rg0'] > 0.7
    assert s.edges_stable > 100
    assert not hasattr(nw, 'ever_fsw')
    return sim, s


# ---------- 2. SWNetwork only via sti.Sim ----------

def test_sw_only():
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sti.SWNetwork()],
        demographics=DEMOG,
        dur=DUR,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    nw = sim.networks.swnetwork
    s = network_stats(nw, sim.people)
    print_summary('SWNetwork only', s)

    # Observed @ N_AGENTS=1000, DUR=10, SEED=1: ever_fsw≈33, ever_client≈81,
    # current_client≈9, lifetime_sw_partners_mean≈1.9.
    # Defaults: fsw_shares=0.05, client_shares=0.12 → for ~500 females and
    # ~500 males alive over the run, expect roughly those raw counts.
    assert 15 < s.ever_fsw_count < 70
    assert 40 < s.ever_client_count < 130
    assert 0 <= s.current_fsw_count <= s.ever_fsw_count
    assert 0 <= s.current_client_count <= s.ever_client_count
    assert s.lifetime_sw_partners_mean > 0.5
    return sim, s


# ---------- 3. MFNetwork + SWNetwork (modular) ----------

def test_mf_plus_sw_modular():
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sti.MFNetwork(), sti.SWNetwork()],
        demographics=DEMOG,
        dur=DUR,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    mf, sw = sim.networks.mfnetwork, sim.networks.swnetwork
    sm = network_stats(mf, sim.people)
    sn = network_stats(sw, sim.people)
    print_summary('MF (modular)', sm)
    print_summary('SW (modular)', sn)

    # Modular MF should match standalone MF closely (same seed, same module);
    # SW similarly. Drift is allowed only from agents added/removed via demog.
    assert 1.0 < sm.lifetime_partners_f_mean < 2.5
    assert sm.edges_stable > 100
    assert 15 < sn.ever_fsw_count < 70
    assert 'fsw' not in [s_.name for s_ in mf.state_list], 'SW state must live on SWNetwork only'
    return sim, sm, sn


# ---------- 4. StructuredSexual (standard MF + SW combined network) ----------

def test_structured_sexual_bundle():
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sti.StructuredSexual()],
        demographics=DEMOG,
        dur=DUR,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    nw = sim.networks.structuredsexual
    s = network_stats(nw, sim.people)
    print_summary('StructuredSexual', s)

    # Observed @ N_AGENTS=1000, DUR=10, SEED=1: f_mean≈1.47, ever_fsw≈21,
    # ever_client≈58, lifetime_sw_partners_mean≈0.65. Slightly lower than the
    # modular case because the bundled network shares an edge list and FSW/MF
    # matching pools draw from the same agents.
    assert 1.0 < s.lifetime_partners_f_mean < 2.5
    assert 0.3 < s.lifetime_stable_partners_mean < 1.0
    assert 10 < s.ever_fsw_count < 60
    assert 30 < s.ever_client_count < 100
    assert s.lifetime_sw_partners_mean > 0.2
    assert set(nw.edge_types) == {'stable', 'casual', 'onetime', 'sw'}
    assert (nw.fsw & ~nw.ever_fsw).sum() == 0
    return sim, s


# ---------- 5. Implicit SW via behavioral thresholds (MF only) ----------

def test_implicit_sw_via_thresholds():
    """No SW network — define 'sex workers' as MF agents whose lifetime
    partner count exceeds a behavioral threshold."""
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sti.MFNetwork(prop_f2=0.05, f2_conc=0.5)],
        demographics=DEMOG,
        dur=DUR,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    nw = sim.networks.mfnetwork
    PARTNER_THRESHOLD = 5
    high_activity_f = (nw.lifetime_partners.values >= PARTNER_THRESHOLD) & sim.people.female.values
    print(f'\n--- Implicit-SW (≥{PARTNER_THRESHOLD} lifetime partners) ---')
    n_high = int(high_activity_f.sum())
    n_female = int(sim.people.female.sum())
    print(f'  high-activity females: {n_high} / {n_female}')
    # Observed @ N_AGENTS=1000, DUR=10, SEED=1, prop_f2=0.05, f2_conc=0.5:
    # ~70 / ~590 — i.e. roughly the bumped high-risk share converted into a
    # behavioral threshold. Should be a non-trivial minority.
    assert 30 < n_high < 200
    assert n_high < 0.3 * n_female
    return sim, high_activity_f


# ---------- 6. SW window dynamics (entry / exit) ----------

def test_sw_window_dynamics():
    """Tight age window — current FSW count should sit below ever_fsw."""
    sw = sti.SWNetwork(
        age_sw_start=ss.normal(loc=22, scale=2),
        dur_sw=ss.lognorm_ex(mean=4, std=1),
    )
    sim = sti.Sim(
        n_agents=N_AGENTS,
        diseases=['hiv'],
        networks=[sw],
        demographics=DEMOG,
        dur=20,
        rand_seed=SEED,
        verbose=0,
    )
    sim.run()
    nw = sim.networks.swnetwork
    s = network_stats(nw, sim.people)
    print_summary('Windowed SW (20-yr run)', s)

    # Observed @ N_AGENTS=1000, DUR=20, SEED=1, mean dur_sw=4y:
    # ever_fsw≈35, current_fsw≈0–3 (most have aged out of their ~4y window).
    assert 15 < s.ever_fsw_count < 80
    assert s.current_fsw_count < s.ever_fsw_count
    # Window narrows the active pool — ratio should be small over a 20y run
    # with a 4y window
    if s.ever_fsw_count > 10:
        ratio = s.current_fsw_count / s.ever_fsw_count
        print(f'  current/ever fsw ratio: {ratio:.2f}')
        assert ratio < 0.5

    # Clamp invariant: age_sw_start >= debut
    ever_uids = nw.ever_fsw.uids
    if len(ever_uids):
        assert (nw.age_sw_start[ever_uids] >= nw.debut[ever_uids]).all()
    return sim, s


# ---------- 7. sti.Sim construction patterns (mirrors test_sim_creation) ----------

def test_sim_construction_patterns():
    """All three sti.Sim creation styles continue to work with the new classes."""

    # Style A — pass a network instance in a list (sti.Sim treats a bare
    # instance as "no networks supplied" and falls back to defaults)
    sim_a = sti.Sim(
        n_agents=300, diseases=['hiv'],
        networks=[sti.StructuredSexual(prop_f0=0.7)],
        demographics=DEMOG, dur=3, rand_seed=SEED, verbose=0,
    )
    sim_a.init()
    assert sim_a.networks.structuredsexual.pars.prop_f0 == 0.7

    # Style B — nw_pars routed by sti.Sim into the default StructuredSexual
    sim_b = sti.Sim(
        pars=dict(n_agents=300, dur=3),
        nw_pars=dict(prop_f0=0.6, fsw_shares=0.1),
        diseases=['hiv'],
        rand_seed=SEED, verbose=0,
    )
    sim_b.init()
    assert sim_b.networks.structuredsexual.pars.prop_f0 == 0.6
    assert float(sim_b.networks.structuredsexual.pars.fsw_shares.pars.p) == 0.1

    # Style C — flat kwargs (mixes sim, network, sti pars)
    sim_c = sti.Sim(
        n_agents=300, dur=3, beta_m2f=0.05, prop_f0=0.5, fsw_shares=0.07,
        diseases=['hiv'], rand_seed=SEED, verbose=0,
    )
    sim_c.init()
    assert sim_c.networks.structuredsexual.pars.prop_f0 == 0.5
    assert sim_c.diseases.hiv.pars.beta_m2f == 0.05

    print('\nsti.Sim construction patterns OK')
    return sim_a, sim_b, sim_c


# ---------- 8. hivsim.Sim defaults ----------

def test_hivsim_sim_defaults():
    sim = hivsim.Sim(n_agents=500, dur=5, rand_seed=SEED, verbose=0)
    sim.run()
    assert 'structuredsexual' in sim.networks
    nw = sim.networks.structuredsexual
    s = network_stats(nw, sim.people)
    print_summary('hivsim.Sim defaults', s)
    # Observed @ n_agents=500, dur=5, SEED=1: f_mean≈1.27, m_mean≈1.27,
    # ever_fsw≈7, ever_client≈28
    assert 0.5 < s.lifetime_partners_f_mean < 2.5
    assert 0.5 < s.lifetime_partners_m_mean < 2.5
    assert 3 < s.ever_fsw_count < 25
    assert 10 < s.ever_client_count < 60
    return sim, s


# ---------- 9. hivsim.demo('simple') uses MFNetwork ----------

def test_hivsim_demo_simple():
    sim = hivsim.demo('simple', run=False, n_agents=500, dur=5, rand_seed=SEED, verbose=0)
    sim.run()
    assert 'mfnetwork' in sim.networks, "simple demo should use MFNetwork"
    assert 'structuredsexual' not in sim.networks
    nw = sim.networks.mfnetwork
    s = network_stats(nw, sim.people)
    print_summary("hivsim.demo('simple') — MFNetwork", s)
    # Observed @ n_agents=500, dur=5, SEED=1: f_mean≈1.29, m_mean≈1.32,
    # mean_dur_stable≈1080, mean_dur_casual≈12.
    assert 0.5 < s.lifetime_partners_f_mean < 2.5
    assert 0.5 < s.lifetime_partners_m_mean < 2.5
    assert 600 < s.mean_dur_stable < 1400
    assert 5 < s.mean_dur_casual < 25
    return sim, s


# ---------- 10. hivsim.demo('zimbabwe') uses StructuredSexual ----------

def test_hivsim_demo_zimbabwe():
    sim = hivsim.demo('zimbabwe', run=False, n_agents=500, stop=1995, rand_seed=SEED, verbose=0)
    sim.run()
    assert 'structuredsexual' in sim.networks, "zimbabwe demo should use StructuredSexual"
    nw = sim.networks.structuredsexual
    s = network_stats(nw, sim.people)
    print_summary("hivsim.demo('zimbabwe') — StructuredSexual", s)
    # Observed @ n_agents=500, stop=1995, SEED=1: f_mean≈0.91, ever_fsw≈13,
    # ever_client≈32. (lifetime_sw_partners_mean is noisy at small N — skip.)
    assert 0.5 < s.lifetime_partners_f_mean < 2.0
    assert 5 < s.ever_fsw_count < 30
    assert 15 < s.ever_client_count < 60
    return sim, s


# ---------- Run-as-script ----------

if __name__ == '__main__':
    out = sc.objdict()
    out.mf_only = test_mf_only()
    out.sw_only = test_sw_only()
    out.mf_plus_sw = test_mf_plus_sw_modular()
    out.bundle = test_structured_sexual_bundle()
    out.implicit = test_implicit_sw_via_thresholds()
    out.window = test_sw_window_dynamics()
    out.creation = test_sim_construction_patterns()
    out.hivsim_default = test_hivsim_sim_defaults()
    out.hivsim_simple = test_hivsim_demo_simple()
    out.hivsim_zimbabwe = test_hivsim_demo_zimbabwe()
    print('\nAll devtests passed.')
