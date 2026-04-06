"""
STI scientific validation: verify that disease parameters affect epi dynamics
in the expected directions and that disease-specific features work correctly.

Covers all non-HIV STIs (syphilis, gonorrhea, chlamydia, trichomoniasis, BV, GUD)
and HIV epi parameter sensitivity. HIV-specific scientific validation (CD4, ART
effects, natural history) lives in test_hiv.py.
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

n_agents = 2000
debug = False


def test_syph_epi():
    """ Test that syphilis epi parameters affect dynamics in the expected direction """
    sc.heading('Test epi dynamics of syphilis')

    base_pars = dict(n_agents=n_agents)

    # Define the parameters to vary
    par_effects = dict(
        rel_trans_primary=[0.5, 1.5],  # Relative transmissibility during primary infection
        init_prev=[0.1, 0.9],
        beta_m2f=[0.2, 0.5]  # Beta for male to female transmission; opposite direction uses half this value
    )

    # Loop over each of the above parameters and make sure they affect the epi dynamics in the expected ways
    for par, par_val in par_effects.items():
        lo = par_val[0]
        hi = par_val[1]

        # Make baseline pars
        pars0 = sc.dcp(base_pars)
        pars1 = sc.dcp(base_pars)

        simpardict_lo = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: lo})
        simpardict_hi = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: hi})

        pars0['diseases'] = sti.Syphilis(**simpardict_lo)
        pars1['diseases'] = sti.Syphilis(**simpardict_hi)

        # Run the simulations and pull out the results
        s0 = sti.Sim(pars0, label=f'{par} {par_val[0]}')
        s1 = sti.Sim(pars1, label=f'{par} {par_val[1]}')
        ss.parallel(s0, s1, debug=debug)

        # Check results
        ind = 1 if par == 'init_prev' else -1
        v0 = s0.results.syph.cum_infections[ind]
        v1 = s1.results.syph.cum_infections[ind]

        print(f'Checking with varying {par:10s} ... ', end='')
        assert v0 <= v1, f'Expected infections to be lower with {par}={lo} than with {par}={hi}, but {v0} > {v1})'
        print(f'✓ ({v0} <= {v1})')

    return s0, s1


def test_hiv_epi():
    """ Test that HIV epi parameters affect dynamics in the expected direction """
    sc.heading('Test epi dynamics of hiv')

    base_pars = dict(n_agents=n_agents, diseases='hiv', beta_m2f=0.05, init_prev=0.05)

    # Define the parameters to vary
    par_effects = dict(
        beta_m2f=[0.01, 0.2],
    )

    # Loop over each of the above parameters and make sure they affect the epi dynamics in the expected ways
    for par, par_val in par_effects.items():
        lo = par_val[0]
        hi = par_val[1]

        # Make baseline pars
        pars0 = sc.dcp(base_pars)
        pars1 = sc.dcp(base_pars)

        pars0[par] = lo
        pars1[par] = hi

        # Run the simulations and pull out the results
        s0 = sti.Sim(pars0, label=f'{par} {par_val[0]}')
        s1 = sti.Sim(pars1, label=f'{par} {par_val[1]}')
        ss.parallel(s0, s1, debug=debug)

        # Check results
        ind = 1 if par == 'init_prev' else -1
        v0 = s0.results.hiv.cum_infections[ind]
        v1 = s1.results.hiv.cum_infections[ind]

        print(f'Checking with varying {par:10s} ... ', end='')
        assert v0 <= v1, f'Expected infections to be lower with {par}={lo} than with {par}={hi}, but {v0} > {v1})'
        print(f'✓ ({v0} <= {v1})')

    return s0, s1


def test_stis(which='discharging', n_agents=5e3, start=2010, stop=2020):
    """ Test running grouped STIs: discharging (NG, CT, TV, BV) and ulcerative (syphilis, GUD) """
    sc.heading('Test STI sim')

    if which == 'discharging':
        sti_pars = dict(
            ng=dict(beta_m2f=0.06, init_prev=0.02),
            ct=dict(beta_m2f=0.06, init_prev=0.05),
            tv=dict(beta_m2f=0.1, init_prev=0.1),
            bv=dict(),
        )
    elif which == 'ulcerative':
        sti_pars = dict(
            syph=dict(beta_m2f=0.1, init_prev=0.01),
            gud=dict(prevalence=0.05),
        )

    pregnancy = ss.Pregnancy(fertility_rate=10)
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


def test_bv(include_hiv=False, n_agents=500):
    """ Test BV dynamics with optional HIV co-infection and a custom menstrual hygiene intervention """

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

    bv = sti.BV(init_prev=0.025)
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
    s1 = ss.Sim(**sim_args, interventions=intvs + [menstrual_hygiene(start=2020, new_val=0.1)])
    ss.parallel(s0, s1, debug=debug)
    return s0, s1


if __name__ == '__main__':
    sc.options(interactive=False)
    s1, s2 = test_syph_epi()
    s3, s4 = test_hiv_epi()
    test_stis(which='discharging')
    test_stis(which='ulcerative')
    test_bv()
    test_bv(include_hiv=True)
