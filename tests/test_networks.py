import stisim as sti
import starsim as ss
import numpy as np
import sciris as sc

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


def test_debut_age():
    """
    Test the debut age in the structured sexual network.
    """

    # Create a structured sexual network with default parameters
    network = sti.StructuredSexual()
    late_debut_network = sti.StructuredSexual(pars={'debut_pars_f': [25, 3], 'debut_pars_m': [26, 3]})
    analyzer = sti.DebutAge()

    s1 = sti.Sim(networks=[network], analyzers=[analyzer])
    s2 = sti.Sim(networks=[late_debut_network], analyzers=[analyzer])

    # Run the simulation
    ss.parallel(s1, s2)

    # all values in the debut age analyzer prop_active_f and prop_active_m should be greater in s1 than in s2
    assert np.all(s1.analyzers.debutage.prop_active_f[0] >= s2.analyzers.debutage.prop_active_f[0]), "Proportion of females active should be higher in default network than in late debut network at any given age"
    assert np.all(s1.analyzers.debutage.prop_active_m[0] >= s2.analyzers.debutage.prop_active_m[0]), "Proportion of males active should be higher in default network than in late debut network at any given age"


if __name__ == '__main__':
    test_network_degrees()
    test_pair_formation()
    test_relationship_duration()
    test_debut_age()
