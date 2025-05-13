"""
Simple sim tests
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti
import hivsim as hs
import pandas as pd
import numpy as np
import pylab as pl

 
def test_hiv_sim(n_agents=500):
    sc.heading('Test simplest possible HIV sim ')

    hiv = hs.HIV(
        beta={'structuredsexual': [0.05, 0.25], 'maternal': [0.05, 0.]},
        init_prev=0.05,
    )
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.FastStructuredSexual()
    maternal = ss.MaternalNet()
    testing = hs.HIVTest(test_prob_data=0.2, start=2000)
    art = hiv.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
    vmmc = hiv.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
    sim = ss.Sim(
        dt=1/12,
        start=1990,
        dur=40,
        n_agents=n_agents,
        diseases=hiv,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=[testing, art, vmmc]
    )
    sim.run(verbose=1/12)

    return sim


def test_msm_hiv(n_agents=500):
    hiv = hs.HIV(beta={'msm': [0.1, 0.1]}, init_prev=0.05)
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    msm = sti.AgeMatchedMSM()
    sim = ss.Sim(
        dt=1/12,
        start=1990,
        dur=10,
        n_agents=n_agents,
        diseases=hiv,
        networks=msm,
        demographics=[pregnancy, death],
    )
    sim.run(verbose=1/12)

    return sim



def test_bv(include_hiv=False, n_agents=500, start=2015, n_years=10):

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
        # p_poor_menstrual_hygiene=0.3,  # Proportion with poor menstrual hygiene
        init_prev=0.025,
    )
    nets = []
    dem = []
    intvs = []
    dis = [bv]
    con = []

    if include_hiv:
        hiv = hiv.HIV(
            beta_m2f=0.5,
            beta_m2c=0.1,
            init_prev=0.1,
        )

        nets += [sti.FastStructuredSexual(store_register=True), ss.MaternalNet()]

        pregnancy = ss.Pregnancy(fertility_rate=10)
        death = ss.Deaths(death_rate=10)
        dem += [pregnancy, death]

        testing = hs.HIVTest(test_prob_data=0.2, start=2000)
        art = hs.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
        vmmc = hs.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
        intvs += [testing, art, vmmc]
        dis += [hiv]
        con += [sti.hiv_bv(hiv_module=hiv, bv_module=bv)]

    # Make sim
    sim_args = dict(unit='year', dt=1/12, start=start, dur=n_years, n_agents=n_agents, diseases=dis, networks=nets, demographics=dem, connectors=con)

    s0 = ss.Sim(**sim_args, interventions=intvs)
    s1 = ss.Sim(**sim_args, interventions=intvs + [menstrual_hygiene(start=2020, new_val=0.1)])
    ss.parallel(s0, s1)
    return [s0, s1]


def test_stis(which='discharging', n_agents=5e3, start=2010, stop=2020):
    sc.heading('Test STI sim')

    if which == 'discharging':
        ng = sti.Gonorrhea(beta_m2f=0.06, init_prev=0.02)
        ct = sti.Chlamydia(beta_m2f=0.06, init_prev=0.05)
        tv = sti.Trichomoniasis(beta_m2f=0.1, init_prev=0.1)
        stis = [ng, ct, tv]
    elif which == 'ulcerative':
        sy = sti.Syphilis(beta_m2f=0.1, init_prev=0.01)
        gud = sti.GUDPlaceholder(prevalence=0.05)
        stis = [sy, gud]

    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.FastStructuredSexual()

    sim = ss.Sim(
        dt=1/12,
        n_agents=n_agents,
        start=start,
        stop=stop,
        diseases=stis,
        networks=sexual,
        demographics=[pregnancy, death],
    )

    sim.run(verbose=1/12)

    return sim


if __name__ == '__main__':

    do_plot = True

    s0 = test_hiv_sim()
    s1 = test_msm_hiv()
    s2 = test_stis(which='discharging')

    if do_plot:
        s1.plot("ng")
        pl.show()

    sims = test_bv(include_hiv=True)
    if do_plot:
        import pylab as pl
        r0 = sims[0].results.bv.prevalence
        r1 = sims[1].results.bv.prevalence
        t = sims[0].results.bv.timevec
        pl.figure()
        pl.plot(t, r0, label='Baseline')
        pl.plot(t, r1, label='Improved menstrual hygiene')
        # pl.axvline(x=2020, color='k', ls='--')
        pl.title('BV prevalence')
        pl.legend()
        pl.show()
