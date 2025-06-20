"""
Test impact of including the HIV-Syphilis connector

One simulation is run with HIV and syphilis, but no connector. A second simulation
is run with the addition of the connector. Adding the connector should increase
transmission of HIV overall due to syphilis coinfection. The prevalance should
therefore be higher when the connector is included.
"""

import sciris as sc
import starsim as ss
import stisim as sti
import pylab as pl
import numpy as np

from stisim import coinfection_stats

np.seterr(all='raise')

def ret(*values):
    return values if __name__ == '__main__' else None

def test_hiv_syph():

    # Make diseases
    hiv = sti.HIV(init_prev=0.1, beta={'structuredsexual': [0.01, 0.01]})
    syphilis = sti.SyphilisPlaceholder(prevalence=0.9)
    coinfection_analyzer = coinfection_stats(disease1=syphilis, disease2='hiv', disease1_infected_state_name='active')

    pars = dict(
        start=2000,
        stop=2020,
        dt=1/12,
        verbose=1/12,
        n_agents=1000,
        analyzers=[coinfection_analyzer],
        networks=sti.StructuredSexual(),
        diseases=[hiv, syphilis]
    )
    s0 = ss.Sim(pars).run()
    r0 = s0.results.hiv.cum_infections[-1]

    pars['connectors'] = sti.hiv_syph(hiv, syphilis, rel_sus_hiv_syph=2, rel_trans_hiv_syph=2)
    s1 = ss.Sim(pars).run()
    r1 = s1.results.hiv.cum_infections[-1]

    assert r0 <= r1, f'The hiv-syph connector should increase HIV infections, but {r1}<{r0}'
    print(f'✓ ({r0} <= {r1})')

    return ret(s0, s1)


# %% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    s0, s1 = test_hiv_syph()

    pl.plot(s0.results.hiv.timevec, s0.results.hiv.prevalence, label='No connector')
    pl.plot(s1.results.hiv.timevec, s1.results.hiv.prevalence, label='With connector')
    pl.ylim(0, 1)
    pl.legend()
    pl.xlabel('Year')
    pl.ylabel('Prevalence')
    pl.show()

    sc.toc(T)
    print('Done.')
