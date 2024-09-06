"""
Testing out product mix for gonorrhea
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd
import numpy as np
import pylab as pl
 

class TestingAlgorithm(sti.STITest):
    def __init__(self, pars=None, treatments=None, diseases=None, **kwargs):
        super().__init__(**kwargs)
        self.default_pars(
        )
        self.update_pars(pars, **kwargs)
        self.treatments = treatments
        self.diseases = diseases
        return

    def apply(self, sim, uids=None):
        """ What needs to happen? """
        # Figure out who gets what product

        return


def test_gon(n_agents=500, dt=1/12, start=2000, n_years=40):

    ng = sti.Gonorrhea(
        beta_m2f=0.038,
        beta_f2m=0.019,
        init_prev=0.01,
    )

    bv = sti.DischargingSTI(
        beta_m2f=0.038,
        beta_f2m=0.019,
        init_prev=0.01,
    )

    # Care seekers interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic
        bv_care = sim.diseases.bv.symptomatic
        return (bv_care | ng_care).uids

    # Okayyyyy
    df = sc.dataframe.read_csv('test_data/ng_dx_product_mix.csv')
    dxprods = dict()
    for name in df.name.unique():
        dxprods[name] = sti.STIDx(df[df.name == name])
    return dxprods

    pm = pd.read_csv('test_data/ng_dx_product_mix.csv')
    product_mix = sti.ProductMix(
        df=pm,
        products=None
    )



    ng_tx = sti.GonorrheaTreatment()

    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual()
    sim = ss.Sim(
        dt=dt,
        start=start,
        n_years=n_years,
        n_agents=n_agents,
        diseases=[ng, bv],
        networks=sexual,
        demographics=[pregnancy, death],
        interventions=[syndromic, ng_tx, tv_tx, ct_tx, vd_tx],
        analyzers=[],
    )
    sim.run(verbose=0.1)

    return sim


if __name__ == '__main__':
    sim = test_gon(n_agents=5e3, dt=1/12, n_years=20)
    sim.plot('gonorrheatreatment')
    pl.show()
