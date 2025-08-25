"""
Simple sim tests
"""
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd
import numpy as np

debug = False  # Run in serial


def test_hiv_sim(n_agents=500):
    sc.heading('Test simplest possible HIV sim ')

    hiv = sti.HIV(
        beta_m2f=0.05,
        beta_m2c=0.1,
        init_prev=0.05,
    )
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual(recall_prior=True)
    prior = sti.PriorPartners()
    maternal = ss.MaternalNet()
    testing = sti.HIVTest(test_prob_data=0.2, start=2000)
    art = sti.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
    vmmc = sti.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
    sim = sti.Sim(
        start=1990,
        dur=40,
        n_agents=n_agents,
        diseases=hiv,
        networks=[sexual, prior, maternal],
        demographics=[pregnancy, death],
        interventions=[testing, art, vmmc]
    )
    sim.run(verbose=1/12)

    return sim


def test_msm_hiv(n_agents=500):
    hiv = sti.HIV(beta_m2m=0.1, init_prev=0.05)
    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    msm = sti.AgeMatchedMSM()
    sim = sti.Sim(
        start=1990,
        dur=10,
        n_agents=n_agents,
        diseases=hiv,
        networks=msm,
        demographics=[pregnancy, death],
    )
    sim.run(verbose=1/12)

    return sim


def test_bv(include_hiv=False, n_agents=500):

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
        hiv = sti.HIV(
            beta_m2f=0.5,
            beta_m2c=0.1,
            init_prev=0.1,
        )

        nets += [sti.StructuredSexual(store_register=True), ss.MaternalNet()]

        pregnancy = ss.Pregnancy(fertility_rate=10)
        death = ss.Deaths(death_rate=10)
        dem += [pregnancy, death]

        testing = sti.HIVTest(test_prob_data=0.2, start=2000)
        art = sti.ART(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_art': np.linspace(0, 0.9, 21)}))
        vmmc = sti.VMMC(coverage_data=pd.DataFrame(index=np.arange(2000, 2021), data={'p_vmmc': np.linspace(0.025, 0.125, 21)}))
        intvs += [testing, art, vmmc]
        dis += [hiv]
        con += [sti.hiv_bv(hiv_module=hiv, bv_module=bv)]

    # Make sim
    sim_args = dict(start=2015, stop=2030, n_agents=n_agents, diseases=dis, networks=nets, demographics=dem, connectors=con)

    s0 = ss.Sim(**sim_args, interventions=intvs)
    s1 = ss.Sim(**sim_args, interventions=intvs + [menstrual_hygiene(start=ss.date(2020), new_val=0.1)])
    ss.parallel(s0, s1, debug=debug)
    return s0, s1


def test_stis(which='discharging', n_agents=5e3, start=2010, stop=2020):
    sc.heading('Test STI sim')

    if which == 'discharging':
        sti_pars = dict(
            ng=dict(beta_m2f=0.06, init_prev=0.02),
            ct=dict(beta_m2f=0.06, init_prev=0.05),
            tv=dict(beta_m2f=0.1, init_prev=0.1),
            bv=dict(),  # Bacterial vaginosis does not have a beta_m2f parameter
        )
    elif which == 'ulcerative':
        sti_pars = dict(
            syph=dict(beta_m2f=0.1, init_prev=0.01),
            gud=dict(prevalence=0.05),  # Placeholder for GUD
        )

    pregnancy = sti.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)

    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        diseases=list(sti_pars.keys()),
        sti_pars=sti_pars,
        demographics=[pregnancy, death],
    )

    sim.run(verbose=1/12)

    return sim


def test_sim_creation():

    start = 2010
    stop = 2020

    pars = dict(
        start=start,
        stop=stop,
    )

    nw_pars = dict(debut=ss.lognorm_ex(20, 5))
    sti_pars = dict(ng=dict(eff_condom=0.6))
    datafolder = './test_data/'

    # Test 1: default networks with custom pars, demographics from location string, and diseases from disease names with custom pars
    sim1 = sti.Sim(
        pars=pars,
        nw_pars=nw_pars,
        demographics='zimbabwe',
        datafolder=datafolder,
        diseases=['ng', 'ct', 'tv', 'bv', 'hiv'],
        sti_pars=sti_pars,
        # connectors=True
    )
    sim1.init()

    assert sim1.diseases.ng.pars.eff_condom == 0.6, "Disease parameter not set correctly"
    assert len(sim1.diseases) == 5, "Incorrect number of diseases initialized"
    # assert len(sim1.connectors) > 0, "No connectors initialized"

    # Test 2: mix of strings and modules
    demographics = [sti.Pregnancy(), ss.Deaths()]  # Replace the default ss.Pregnancy module with the sti one
    networks = sti.StructuredSexual()
    diseases = [sti.Gonorrhea(), 'hiv']

    sim2 = sti.Sim(
        pars=pars,
        networks=networks,
        demographics=demographics,
        diseases=diseases,
        # connectors=True,
    )

    sim2.init()

    assert isinstance(sim2.networks.structuredsexual, sti.StructuredSexual), "Network not initialized correctly"
    assert len(sim2.diseases) == 2, "Incorrect number of diseases initialized"
    # assert len(sim2.connectors) > 0, "No connectors initialized"
    assert len(sim2.demographics) == 2, "Incorrect number of demographics initialized"

    # Test 3: flat pars dict
    pars = dict(
        start=2010,  # Sim par
        beta_m2f=0.05,  # STI parameter applied to all STIs
        prop_f0=0.45,
        location='zimbabwe',
        datafolder='./test_data/',
        diseases=['ng', 'ct', 'tv'],
        ng=dict(eff_condom=0.6),  # Gonorrhea-specific parameter
    )

    sim3 = sti.Sim(**pars)
    sim3.init()

    assert sim3.diseases.ng.pars.beta_m2f == pars['beta_m2f'], "Disease parameter not set correctly"
    assert sim3.diseases.ct.pars.beta_m2f == pars['beta_m2f'], "Disease parameter not set correctly"
    assert sim3.diseases.ng.pars.eff_condom == pars['ng']['eff_condom'], "Disease parameter not set correctly"
    assert sim3.networks.structuredsexual.pars.prop_f0 == pars['prop_f0'], "Network parameter not set correctly"
    assert len(sim3.networks) == 2, "Default networks not added"
    assert len(sim3.diseases) == 3, "Incorrect number of diseases initialized"

    return


def devtest_location():
    """
    Won't currently run on GH actions, but can run locally to check authentication key
    """
    sc.heading('Test location-based sim creation')
    sim1 = sti.Sim(location='zambia', start=1990, stop=2040)
    sim1.run()
    assert len(sim1.demographics) == 3, "Demographics not initialized"

    return


def test_time():
    sim = sti.Sim(start=2010, diseases='hiv')
    sim.run()
    assert sim.pars.start == sim.t.yearvec[0], "Timevec seems incorrect"
    return sim


if __name__ == '__main__':

    do_plot = False

    s0 = test_hiv_sim()
    s1 = test_msm_hiv()
    s2 = test_bv()
    s3 = test_stis(which='discharging')
    test_sim_creation()
    devtest_location()
    s4 = test_time()

