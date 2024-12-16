"""
Simple sim tests
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd
import numpy as np
import pylab as pl

 
def test_hiv_sim(n_agents=500, dt=1):
    sc.heading('Test simplest possible HIV sim ')

    hiv = sti.HIV(
        beta={'structuredsexual': [0.05, 0.25], 'maternal': [0.05, 0.]},
        init_prev=0.05,
    )
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual()
    maternal = ss.MaternalNet()
    testing = sti.HIVTest(test_prob_data=0.2, start=2000)
    art = sti.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
    vmmc = sti.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
    sim = ss.Sim(
        dt=dt,
        start=1990,
        dur=40,
        n_agents=n_agents,
        diseases=hiv,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=[testing, art, vmmc]
    )
    sim.run()

    return sim


def test_sti_sim(n_agents=500, dt=1, start=2000, dur=40):

    bv = sti.DischargingSTI(
        beta_m2f=0.1,
        init_prev=0.025,
    )
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual()
    sim = ss.Sim(
        dt=dt,
        start=start,
        dur=dur,
        n_agents=n_agents,
        diseases=bv,
        networks=sexual,
        demographics=[pregnancy, death],
        analyzers=[],
    )
    sim.run(verbose=0.01)

    return sim


if __name__ == '__main__':
    s0 = test_hiv_sim()
    s1 = test_sti_sim(n_agents=5e3, dt=1/12, dur=20)
    # sim.plot("dischargingsti")
    # pl.show()


