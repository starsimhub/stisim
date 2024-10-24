"""
Test utils
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti

 
def test_sim_results():
    sc.heading('Test aggregating sim results')
    sim = ss.demo(plot=False, verbose=-1, dt=1/4, summary=False)
    df = sti.finalize_results(sim, modules_to_drop=['randomnet'])

    return sim, df


if __name__ == '__main__':
    sim, df = test_sim_results()


