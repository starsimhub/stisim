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


def test_gon(n_agents=500, dt=1, start=2000, n_years=40):

    gon = sti.Gonorrhea(
        beta_m2f=0.038,
        beta_f2m=0.019,
        init_prev=0.03,
    )
    trich = sti.Trichomoniasis(
        beta_m2f=0.05,
        beta_f2m=0.02,
        init_prev=0.05,
    )
    vd = sti.VDPlaceholder(prevalence=0.05)

    # Testing interventions
    def vd_symptoms(sim):
        return (sim.diseases.gonorrhea.symptomatic | sim.diseases.trichomoniasis.symptomatic | sim.diseases.vd.symptomatic).uids

    gon_tx = sti.GonorrheaTreatment()
    syndromic = sti.SyndromicMgmt(
        diseases=['gonorrhea', 'trichomoniasis'],
        test_prob_data=0.2,
        eligibility=vd_symptoms,
        treatment=dict(
            gonorrhea=gon_tx,
        )
    )

    pregnancy = ss.Pregnancy(fertility_rate=10)
    death = ss.Deaths(death_rate=10)
    sexual = sti.StructuredSexual()
    sim = ss.Sim(
        dt=dt,
        start=start,
        n_years=n_years,
        n_agents=n_agents,
        diseases=[gon, trich, vd],
        networks=sexual,
        demographics=[pregnancy, death],
        interventions=[syndromic, gon_tx],
        analyzers=[],
    )
    sim.run(verbose=0.01)

    return sim


if __name__ == '__main__':
    sim = test_gon(n_agents=5e3, dt=1/52, n_years=20)
    sim.plot('gonorrhea')
    pl.show()
