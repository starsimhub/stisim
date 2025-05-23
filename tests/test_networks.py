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


if __name__ == '__main__':
    test_network_degrees()