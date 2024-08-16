"""
Testing out AMR implementations
"""

# Imports
import sciris as sc
import starsim as ss
import stisim as sti
import pandas as pd
import numpy as np
import pylab as pl


def test_gon(n_agents=500, dt=1/12, start=2000, n_years=40):

    ng = sti.Gonorrhea(
        beta_m2f=0.038,
        beta_f2m=0.019,
        init_prev=0.01,
    )
    tv = sti.Trichomoniasis(
        beta_m2f=0.02,
        beta_f2m=0.01,
        init_prev=0.05,
    )
    vd = sti.DischargingSTI(
        beta_m2f=0.1,
        beta_f2m=0.05,
        init_prev=0.025,
    )
    ct = sti.Chlamydia(
        beta_m2f=0.05,
        beta_f2m=0.025,
        init_prev=0.04,
    )

    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.ti)
        tv_care = sim.diseases.tv.symptomatic & (sim.diseases.tv.ti_seeks_care == sim.ti)
        ct_care = sim.diseases.ct.symptomatic & (sim.diseases.ct.ti_seeks_care == sim.ti)
        vd_care = sim.diseases.vd.symptomatic & (sim.diseases.vd.ti_seeks_care == sim.ti)
        return (ng_care | tv_care | ct_care | vd_care).uids

    ng_tx = sti.GonorrheaTreatment()
    tv_tx = sti.STITreatment(disease='tv', name='tv_tx', label='tv_tx')
    ct_tx = sti.STITreatment(disease='ct', name='ct_tx', label='ct_tx')
    vd_tx = sti.STITreatment(disease='vd', name='vd_tx', label='vd_tx')
    syndromic = sti.SyndromicMgmt(
        diseases=[ng, tv, ct, vd],
        eligibility=seeking_care_discharge,
        treatments=[ng_tx, tv_tx, ct_tx, vd_tx],
    )

    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual()
    sim = ss.Sim(
        dt=dt,
        start=start,
        n_years=n_years,
        n_agents=n_agents,
        diseases=[ng, tv, ct, vd],
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
