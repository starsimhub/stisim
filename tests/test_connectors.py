"""
Multi-disease interaction validation: verify that connectors modify
transmission and disease progression as expected.

One simulation is run with HIV and syphilis, but no connector. A second simulation
is run with the addition of the connector. Adding the connector should increase
transmission of HIV overall due to syphilis coinfection. The prevalence should
therefore be higher when the connector is included.
"""

import sciris as sc
import starsim as ss
import stisim as sti
import pylab as pl

from stisim import coinfection_stats

# np.seterr(all='raise')
debug = False


def test_hiv_syph():

    # Make diseases
    def make_args():
        hiv = sti.HIV(init_prev=0.1, beta_m2f=0.05)
        syphilis = sti.SyphilisPlaceholder(prevalence=0.9)
        coinfection_analyzer = coinfection_stats(disease1=syphilis, disease2='hiv', disease1_infected_state_name='active')
        return hiv, syphilis, coinfection_analyzer

    hiv, syphilis, coinfection_analyzer = make_args()
    pars = dict(analyzers=coinfection_analyzer, diseases=[hiv, syphilis], n_agents=500, dur=10)
    s0 = sti.Sim(pars)

    pars['connectors'] = sti.hiv_syph(hiv, syphilis, rel_sus_hiv_syph=20, rel_trans_hiv_syph=20)
    hiv, syphilis, coinfection_analyzer = make_args()
    pars = dict(analyzers=coinfection_analyzer, diseases=[hiv, syphilis], n_agents=500, dur=10)
    s1 = sti.Sim(pars)

    ss.parallel(s0, s1, debug=debug)

    r0 = s0.results.hiv.cum_infections[-1]
    r1 = s1.results.hiv.cum_infections[-1]

    assert r0 <= r1, f'The hiv-syph connector should increase HIV infections, but {r1}<{r0}'
    print(f'✓ hiv-syph connector increased HIV infections ({r0} <= {r1})')

    return s0, s1


def test_single_connector_instance():
    """
    Regression: passing a single connector instance (not a list) via
    ``connectors=`` should not trigger the auto-add machinery — previously
    ``len(connector_instance)`` returned 0 and a duplicate was added,
    causing ndict to raise on init.
    """
    hiv  = sti.HIV(init_prev=0.05, beta_m2f=0.05)
    syph = sti.Syphilis(init_prev=0.05, beta_m2f=0.05)
    sim = sti.Sim(
        diseases=[hiv, syph],
        connectors=sti.hiv_syph(hiv, syph),  # single instance, not a list
        n_agents=200, dur=2, verbose=-1,
    )
    sim.run()
    assert len(sim.pars['connectors']) == 1, \
        f'Expected exactly one connector, got {len(sim.pars["connectors"])}'
    return sim


# %% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    s0, s1 = test_hiv_syph()
    test_single_connector_instance()

    pl.plot(s0.results.hiv.timevec, s0.results.hiv.prevalence, label='No connector')
    pl.plot(s1.results.hiv.timevec, s1.results.hiv.prevalence, label='With connector')
    pl.ylim(0, 1)
    pl.legend()
    pl.xlabel('Year')
    pl.ylabel('Prevalence')
    pl.show()

    sc.toc(T)
    print('Done.')
