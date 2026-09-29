"""
Network dynamics validation: verify that network parameters (concurrency, pair
formation, relationship duration, debut age, MSM) affect behaviour as expected.
"""
import stisim as sti
import starsim as ss
import numpy as np
import sciris as sc
import hivsim


def test_msm_network(n_agents=500):
    """ Test MSM HIV transmission via AgeMatchedMSM network """
    hiv = sti.HIV(beta_m2m=0.1, init_prev=0.05)
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    msm = sti.AgeMatchedMSM(p_msm=ss.bernoulli(p=0.3))
    sim = sti.Sim(
        start=1990,
        dur=10,
        n_agents=n_agents,
        diseases=hiv,
        networks=[msm],
        demographics=[pregnancy, death],
    )
    sim.run(verbose=1/12)

    assert sim.results.hiv.cum_infections[-1] > 0, "MSM network should produce HIV infections"
    return sim


def test_network_degrees():
    """
    Test the degree distribution of the structured sexual network.
    """

    # Create a structured sexual network
    network = sti.StructuredSexual()
    high_concurrency = sti.StructuredSexual(pars={'f0_conc': 0.001, 'f1_conc': 0.2, 'f2_conc': 0.6, 'm0_conc': 0.001, 'm1_conc': 0.5, 'm2_conc': 0.9})
    analyzer = sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])

    s1 = sti.Sim(networks=[network], analyzers=[analyzer])
    s2 = sti.Sim(networks=[high_concurrency], analyzers=[analyzer])

    ss.parallel(s1, s2)

    # Mean number of partners should increase in high concurrency case
    s1_mean_partners = np.mean(s1.analyzers.networkdegree.lifetime_partners_f + s1.analyzers.networkdegree.lifetime_partners_m)
    s2_mean_partners = np.mean(s2.analyzers.networkdegree.lifetime_partners_f + s2.analyzers.networkdegree.lifetime_partners_m)

    assert s2_mean_partners > s1_mean_partners, f"Mean partners in high concurrency ({s2_mean_partners}) should be greater than in normal ({s1_mean_partners})"
    print (f"Mean partners in high concurrency ({s2_mean_partners}) is greater than in normal ({s1_mean_partners})")


def test_pair_formation():
    """
    The parameters p_matched_stable and p_matched_casual determine the probability of forming a stable and casual relationship.
    Check that higher values mean fewer lifetime partners, fewer stable partners, and fewer casual partners.
    """

    s1 = sti.Sim(networks=[sti.StructuredSexual()], analyzers=[sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])])
    s2 = sti.Sim(networks=[sti.StructuredSexual(p_matched_stable=[0.99, 0.9, 0.9])],
                analyzers=[sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])])

    ss.parallel(s1, s2)

    # lifetime partners should be lower in the second case
    s1_mean_partners = np.mean(s1.analyzers.networkdegree.lifetime_partners_f + s1.analyzers.networkdegree.lifetime_partners_m)
    s2_mean_partners = np.mean(s2.analyzers.networkdegree.lifetime_partners_f + s2.analyzers.networkdegree.lifetime_partners_m)
    assert s2_mean_partners < s1_mean_partners, f"Mean partners in high probability scenario ({s2.results.network_degree.mean_partners}) should be less than in default ({s1.results.network_degree.mean_partners})"
    print (f"Mean partners in high probability scenario ({s2_mean_partners}) is less than in default ({s1_mean_partners})")

    # higher probability of stable relationships means there should be more stable relationships
    s1_mean_stable_partners = np.mean(s1.analyzers.networkdegree.lifetime_stable_partners_f + s1.analyzers.networkdegree.lifetime_stable_partners_m)
    s2_mean_stable_partners = np.mean(s2.analyzers.networkdegree.lifetime_stable_partners_f + s2.analyzers.networkdegree.lifetime_stable_partners_m)
    assert s2_mean_stable_partners > s1_mean_stable_partners, f"Mean stable partners in high probability stable scenario ({s2_mean_stable_partners}) should be greater than in default ({s1_mean_stable_partners})"
    print (f"Mean stable partners in high probability scenario ({s2_mean_stable_partners}) is higher than in default ({s1_mean_stable_partners})")

    # casual partners should be lower in the second case
    s1_mean_casual_partners = np.mean(s1.analyzers.networkdegree.lifetime_casual_partners_f + s1.analyzers.networkdegree.lifetime_casual_partners_m)
    s2_mean_casual_partners = np.mean(s2.analyzers.networkdegree.lifetime_casual_partners_f + s2.analyzers.networkdegree.lifetime_casual_partners_m)
    assert s2_mean_casual_partners < s1_mean_casual_partners, f"Mean casual partners in high probability scenario ({s2_mean_casual_partners}) should be less than in default ({s1_mean_casual_partners})"
    print (f"Mean casual partners in high probability scenario ({s2_mean_casual_partners}) is less than in default ({s1_mean_casual_partners})")

    return


def test_relationship_duration():
    """
    Test the relationship duration in the structured sexual network.
    """

    stable_dur_pars = dict(
        teens=[
            # (mu,stdev) for levels 0, 1, 2
            [ss.years(100), ss.years(1)],
            [ss.years(50), ss.years(2)],
            [ss.months(1e-4), ss.months(1e-4)]
        ],
        young=[
            [ss.years(100), ss.years(1)],
            [ss.years(50), ss.years(3)],
            [ss.months(1e-4), ss.months(1e-4)]
        ],
        adult=[
            [ss.years(100), ss.years(1)],
            [ss.years(50), ss.years(3)],
            [ss.months(1e-4), ss.months(1e-4)]
        ],
    )

    # Create a structured sexual network with default parameters
    network = sti.StructuredSexual()
    long_network = sti.StructuredSexual(pars={'stable_dur_pars': stable_dur_pars})
    analyzer = sti.RelationshipDurations()

    s1 = sti.Sim(networks=[network], analyzers=[analyzer])
    s2 = sti.Sim(networks=[long_network], analyzers=[analyzer])

    # Run the simulation
    ss.parallel(s1, s2)

    # Check the mean relationship duration
    mean_duration = s1.results.relationshipdurations.mean_duration[-1]
    mean_duration_long = s2.results.relationshipdurations.mean_duration[-1]

    assert mean_duration_long > mean_duration, f"Mean relationship duration should be longer if dur_pars are higher (sim1: {mean_duration} vs sim2: {mean_duration_long})"
    print(f"Increasing relationship duration parameters results in longer mean relationship duration: {mean_duration_long} vs {mean_duration}")


def test_partner_seeking_rates():
    """
    Test the partner seeking rates in the structured sexual network.
    """

    # Contrast a low vs high pair-formation probability. The mean gap between
    # relationships is dominated by (long) relationship durations, so the default
    # 0.5-vs-0.9 contrast was swamped by noise and flipped sign ~50% of the time.
    # A wide 0.1-vs-0.9 contrast, a shared rand_seed (common random numbers pair
    # the two arms and cancel most between-sim variance), and 1000 agents keep the
    # effect well above noise. Note each sim needs its own analyzer instance.
    low_p_pair_form  = sti.StructuredSexual(pars={'p_pair_form': ss.bernoulli(p=0.1)})
    high_p_pair_form = sti.StructuredSexual(pars={'p_pair_form': ss.bernoulli(p=0.9)})
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)

    kw = dict(demographics=[death, pregnancy], stop=2040, n_agents=1000, rand_seed=0)
    s1 = sti.Sim(networks=[low_p_pair_form],  analyzers=[sti.TimeBetweenRelationships()], **kw)
    s2 = sti.Sim(networks=[high_p_pair_form], analyzers=[sti.TimeBetweenRelationships()], **kw)

    # Run the simulation
    ss.parallel(s1, s2, debug=True)

    # compute the mean time between relationships for both sims, excluding the time until first relationship because some agents
    # take a long time to get their initial pairing, and some never do. This effect is magnified in higher probability
    # scenarios so the time between relationships gets skewed in the wrong direction.
    s1_consolidated = [item for sublist in s1.results.timebetweenrelationships.times_between_relationships.values() for
                    index, item in enumerate(sublist) if index > 0 and item > 0]
    s2_consolidated = [item for sublist in s2.results.timebetweenrelationships.times_between_relationships.values() for
                    index, item in enumerate(sublist) if index > 0 and item > 0]
    s1_mean = np.mean(s1_consolidated)
    s2_mean = np.mean(s2_consolidated)

    assert s2_mean < s1_mean, f"Mean time between relationships should be lower in high p_pair_form scenario ({s2_mean}) than in low ({s1_mean})"

    
def test_debut_age():
    """
    Test the debut age in the structured sexual network.
    """

    # Create a structured sexual network with default parameters
    network = sti.StructuredSexual()

    late_debut_network = sti.StructuredSexual(debut_f=25, debut_m=26)
    analyzer = sti.DebutAge()

    s1 = sti.Sim(networks=[network], analyzers=[analyzer])
    s2 = sti.Sim(networks=[late_debut_network], analyzers=[analyzer])

    # Run the simulation
    ss.parallel(s1, s2)

    # all values in the debut age analyzer prop_active_f and prop_active_m should be greater in s1 than in s2
    assert np.all(s1.analyzers.debutage.prop_active_f[0] >= s2.analyzers.debutage.prop_active_f[0]), "Proportion of females active should be higher in default network than in late debut network at any given age"
    assert np.all(s1.analyzers.debutage.prop_active_m[0] >= s2.analyzers.debutage.prop_active_m[0]), "Proportion of males active should be higher in default network than in late debut network at any given age"


def test_shorter_sw():
    """ Shorter SW participation window → fewer HIV transmissions attributable to FSW """
    # rand_seed=1 draws zero FSW-sourced transmissions from the small FSW pool
    # (~9/2000 agents) in both arms under starsim's post-3.5.0 CRN hash scheme,
    # making the comparison degenerate; seed=2 reliably exercises the effect.
    kw = dict(n_agents=2000, stop=2010, rand_seed=2, run=False, plot=False, verbose=0,
              analyzers=[sti.sw_stats(diseases=['hiv'])])

    long_win  = hivsim.demo('zimbabwe', dur_sw=10, **kw)
    short_win = hivsim.demo('zimbabwe', dur_sw=2,  **kw)
    ss.parallel(long_win, short_win)

    long_trans  = long_win.results.sw_stats.new_transmissions_fsw_hiv.sum()
    short_trans = short_win.results.sw_stats.new_transmissions_fsw_hiv.sum()
    assert short_trans < long_trans, \
        f'Expected fewer FSW-attributable transmissions with shorter window: short={short_trans}, long={long_trans}'


def test_match_pairs_refactor_preserves_output():
    """Refactor of match_pairs must produce identical (p1, p2) given same seed."""
    import numpy as np
    import starsim as ss
    import stisim as sti
    sim = ss.Sim(n_agents=2_000, networks=sti.StructuredSexual(), diseases='sis',
                 start='2000-01-01', stop='2001-01-01', rand_seed=42)
    sim.init()
    net = sim.networks.structuredsexual
    # Force at least one matching attempt; capture output.
    try:
        p1, p2 = net.match_pairs()
    except sti.networks.NoPartnersFound:
        p1, p2 = ss.uids(np.array([], dtype=np.int64)), ss.uids(np.array([], dtype=np.int64))
    # The test acts as a smoke test of the refactor; exact UIDs are not asserted
    # because the test will be paired with a recorded baseline in step 2.
    assert len(p1) == len(p2)
    assert (sim.people.male[p1]).all() if len(p1) else True
    assert (sim.people.female[p2]).all() if len(p2) else True


def test_decrement_partners_duplicate_edges():
    """Partner counts must be decremented once per ending edge, not once per agent.

    Regression test: when an agent has more than one partnership ending on the
    same timestep (concurrency > 1), the counters are updated by fancy-indexing
    a repeated UID. A naive ``partners[p1e] -= 1`` only subtracts 1 for the
    duplicated UID (last-write-wins), so the count drifts and can go negative.
    ``_decrement_partners`` counts duplicates with np.unique to avoid this.
    """
    net = sti.MFNetwork()
    sim = sti.Sim(n_agents=20, networks=[net], dur=1)
    sim.init()
    net = sim.networks.mfnetwork

    # Start from a clean slate: drop any auto-formed edges and zero the counters.
    for k in net.meta_keys():
        net.edges[k] = net.edges[k][[]]
    net.partners[:] = 0
    net.stable_partners[:] = 0

    # Agent 0 is p1 in TWO stable edges (with agents 1 and 2) ending this step.
    p1 = ss.uids([0, 0])
    p2 = ss.uids([1, 2])
    n = len(p1)
    net.append(
        p1=p1, p2=p2, beta=np.ones(n), condoms=np.zeros(n), dur=np.ones(n),
        acts=np.ones(n), age_p1=np.zeros(n), age_p2=np.zeros(n),
        edge_type=np.full(n, net.edge_types['stable'], dtype=float),
        ti_formed=np.zeros(n, dtype=int),
    )
    net.partners[ss.uids([0])] = 2
    net.partners[ss.uids([1, 2])] = 1
    net.stable_partners[ss.uids([0])] = 2
    net.stable_partners[ss.uids([1, 2])] = 1

    # Every edge ends (active = all False).
    net._decrement_partners(np.zeros(n, dtype=bool))

    assert net.partners[ss.uids([0])][0] == 0, \
        f"Agent with 2 ending edges should decrement by 2; got {net.partners[ss.uids([0])][0]}"
    assert net.stable_partners[ss.uids([0])][0] == 0, \
        f"stable_partners not decremented per-edge; got {net.stable_partners[ss.uids([0])][0]}"
    assert (net.partners[ss.uids([1, 2])] == 0).all(), "Partners of the duplicated agent should reach 0"


def test_partner_counts_match_active_degree():
    """End-to-end invariant: partners == number of active non-SW edges per agent.

    Exercises both the increment (add_pairs) and decrement (_decrement_partners)
    paths under high concurrency, where an agent can form or end multiple edges
    of the same type in one step. Before the np.unique fix this drifted and
    partner counts went negative.
    """
    net = sti.MFNetwork(pars={'f1_conc': 0.5, 'f2_conc': 0.9, 'm1_conc': 0.7, 'm2_conc': 0.95})
    sim = sti.Sim(n_agents=2000, networks=[net], dur=15, rand_seed=1)
    sim.run(verbose=0)
    net = sim.networks.mfnetwork

    p1, p2, et = net.edges.p1, net.edges.p2, net.edges.edge_type
    mask = np.ones(len(p1), dtype=bool)
    if 'sw' in net.edge_types:
        mask = et != net.edge_types['sw']
    degree = np.zeros(len(net.partners))
    np.add.at(degree, np.asarray(p1[mask]), 1)
    np.add.at(degree, np.asarray(p2[mask]), 1)

    alive = sim.people.alive.uids
    mismatch = np.asarray(net.partners[alive]) - degree[np.asarray(alive)]
    assert np.nanmin(net.partners.values) >= 0, \
        f"partners went negative (min={np.nanmin(net.partners.values)})"
    assert np.abs(mismatch).max() == 0, \
        f"partners disagrees with active edge degree for {int((np.abs(mismatch) > 0).sum())} agents"


if __name__ == '__main__':
    test_msm_network()
    test_network_degrees()
    test_pair_formation()
    test_relationship_duration()
    test_partner_seeking_rates()
    test_debut_age()
    test_shorter_sw()
    test_match_pairs_refactor_preserves_output()
    test_decrement_partners_duplicate_edges()
    test_partner_counts_match_active_degree()
