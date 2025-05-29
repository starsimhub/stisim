import stisim as sti
import starsim as ss


def test_debut_age_analyzer():
    """Test the debut age analyzer."""
    network = sti.StructuredSexual()
    analyzer = sti.DebutAge()

    sim = ss.Sim(networks=[network], analyzers=[analyzer])
    sim.run()
    print("done")



if __name__ == "__main__":
    test_debut_age_analyzer()
    print("Test completed successfully.")