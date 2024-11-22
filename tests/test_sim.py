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


def test_sti_sim(n_agents=500, start=2000, n_years=20):

    bv = sti.BV(
        beta_m2f=0.1,
        init_prev=0.025,
        include_care=False,
    )
    hiv = sti.HIV(
        beta_m2f=0.1,
        beta_m2c=0.03,
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
        unit='year',
        dt=1/12,
        start=start,
        dur=n_years,
        n_agents=n_agents,
        diseases=[bv, hiv],
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=[testing, art, vmmc]
    )
    sim.run(verbose=1/12)

    return sim


def test_bv(n_agents=500, start=2015, n_years=10):

    class menstrual_hygiene(ss.Intervention):
        def __init__(self, pars=None, **kwargs):
            super().__init__()
            self.define_pars(
                start=2000,
                new_val=None,
            )
            self.update_pars(pars, **kwargs)
            return

        def step(self):
            if self.t.now() == self.pars.start:
                self.sim.diseases.bv.pars.p_poor_menstrual_hygiene.set(self.pars.new_val)
            return

    bv = sti.BV(
        beta_m2f=0,  # Don't include sexual transmission
        p_poor_menstrual_hygiene=0.3,  # Proportion with poor menstrual hygiene
        init_prev=0.025,
        include_care=False,
    )
    hiv = sti.HIV(
        beta_m2f=0.1,
        beta_m2c=0.03,
        init_prev=0.05,
    )
    sexual = sti.FastStructuredSexual()
    maternal = ss.MaternalNet()
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    testing = sti.HIVTest(test_prob_data=0.2, start=2000)
    art = sti.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
    vmmc = sti.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
    intvs = [testing, art, vmmc]
    nets = [sexual, maternal]
    dem = [pregnancy, death]
    dis = [bv, hiv]
    sim_args = dict(unit='year', dt=1/12, start=start, dur=n_years, n_agents=n_agents, diseases=dis, networks=nets, demographics=dem, connectors=sti.hiv_bv(hiv_module=hiv, bv_module=bv))

    s0 = ss.Sim(**sim_args, interventions=intvs)
    s1 = ss.Sim(**sim_args, interventions=intvs + [menstrual_hygiene(start=2020, new_val=0.1)])
    ss.parallel(s0, s1)

    return [s0, s1]



if __name__ == '__main__':
    # s0 = test_hiv_sim()
    # s1 = test_sti_sim(n_agents=5e3, n_years=20)
    # s1.plot("bv")
    # pl.show()
    sims = test_bv()

    import pylab as pl
    r0 = sims[0].results.bv.female_adult_prevalence
    r1 = sims[1].results.bv.female_adult_prevalence
    t = sims[0].results.bv.timevec
    pl.figure()
    pl.plot(t, r0, label='Baseline')
    pl.plot(t, r1, label='Improved menstrual hygiene')
    # pl.axvline(x=2020, color='k', ls='--')
    pl.title('BV prevalence')
    pl.legend()
    pl.show()
