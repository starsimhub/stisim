import stisim as sti
import starsim as ss
import numpy as np

def test_network_degrees():
    """
    Test the degree distribution of the structured sexual network.
    """

    # Create a structured sexual network
    network = sti.StructuredSexual()
    high_concurrency = sti.StructuredSexual(pars={'f0_conc': 0.001, 'f1_conc': 0.1, 'f2_conc': 0.2, 'm0_conc': 0.001, 'm1_conc': 0.4, 'm2_conc': 0.9})
    analyzer = sti.NetworkDegree()

    s1 = ss.Sim(networks=[network], analyzers=[analyzer])
    s2 = ss.Sim(networks=[high_concurrency], analyzers=[analyzer])

    ss.parallel(s1, s2)

    # Mean number of partners should increase in high concurrency case
    s1_mean_partners = np.mean(s1.analyzers.networkdegree.lifetime_partners_f + s1.analyzers.networkdegree.lifetime_partners_m)
    s2_mean_partners = np.mean(s2.analyzers.networkdegree.lifetime_partners_f + s2.analyzers.networkdegree.lifetime_partners_m)

    assert s2_mean_partners > s1_mean_partners, f"Mean partners in high concurrency ({s2.results.network_degree.mean_partners}) should be greater than in normal ({s1.results.network_degree.mean_partners})"
    print (f"Mean partners in high concurrency ({s2_mean_partners}) is greater than in normal ({s1_mean_partners})")

def test_pair_formation():
    """
    The parameters p_matched_stable and p_matched_casual determine the probability of forming a stable and casual relationship.
    Check that higher values mean fewer lifetime partners, fewer stable partners, and fewer casual partners.
    """

    s1 = ss.Sim(networks=[sti.StructuredSexual()], analyzers=[sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])])
    s2 = ss.Sim(networks=[sti.StructuredSexual(p_matched_stable=[0.99, 0.9, 0.9], p_mismatched_casual=[0.6, 0.6, 0.6])],
                analyzers=[sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])])

    ss.parallel(s1, s2)

    # lifetime partners should be lower in the second case
    s1_mean_partners = np.mean(s1.analyzers.networkdegree.lifetime_partners_f + s1.analyzers.networkdegree.lifetime_partners_m)
    s2_mean_partners = np.mean(s2.analyzers.networkdegree.lifetime_partners_f + s2.analyzers.networkdegree.lifetime_partners_m)
    assert s2_mean_partners < s1_mean_partners, f"Mean partners in high probability scenario ({s2.results.network_degree.mean_partners}) should be less than in default ({s1.results.network_degree.mean_partners})"
    print (f"Mean partners in high probability scenario ({s2_mean_partners}) is less than in default ({s1_mean_partners})")

    # stable partners should be lower in the second case
    s1_mean_stable_partners = np.mean(s1.analyzers.networkdegree.lifetime_stable_partners_f + s1.analyzers.networkdegree.lifetime_stable_partners_m)
    s2_mean_stable_partners = np.mean(s2.analyzers.networkdegree.lifetime_stable_partners_f + s2.analyzers.networkdegree.lifetime_stable_partners_m)
    assert s2_mean_stable_partners < s1_mean_stable_partners, f"Mean stable partners in high probability scenario ({s2_mean_stable_partners}) should be less than in default ({s1_mean_stable_partners})"
    print (f"Mean stable partners in high probability scenario ({s2_mean_stable_partners}) is less than in default ({s1_mean_stable_partners})")

    # casual partners should be lower in the second case
    s1_mean_casual_partners = np.mean(s1.analyzers.networkdegree.lifetime_casual_partners_f + s1.analyzers.networkdegree.lifetime_casual_partners_m)
    s2_mean_casual_partners = np.mean(s2.analyzers.networkdegree.lifetime_casual_partners_f + s2.analyzers.networkdegree.lifetime_casual_partners_m)
    assert s2_mean_casual_partners < s1_mean_casual_partners, f"Mean casual partners in high probability scenario ({s2_mean_casual_partners}) should be less than in default ({s1_mean_casual_partners})"
    print (f"Mean casual partners in high probability scenario ({s2_mean_casual_partners}) is less than in default ({s1_mean_casual_partners})")

    return







if __name__ == '__main__':
    test_network_degrees()
    test_pair_formation()