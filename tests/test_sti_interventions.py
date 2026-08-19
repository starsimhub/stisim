"""
Test STI testing interventions: STITest/SyphTest eligibility and coverage input formats

Tests cover:
    - STITest.get_testers works with the default eligibility (None -> everyone), not just
      an explicit eligibility function (regression test for GH #537)
    - SyphTest accepts a DataFrame of test probabilities stratified by
      risk group/sex/sex-work status, as documented in the interventions user guide
"""

import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

n_agents = 500


def make_syph_dx():
    df = pd.DataFrame({'state': ['primary', 'secondary', 'latent'], 'p_positive': [0.85, 0.98, 0.7]})
    return sti.SyphDx(df)


def test_stitest_default_eligibility():
    """ STITest/SyphTest must not crash when eligibility is left at its default (None) """
    sc.heading('Testing STITest with default eligibility...')

    test = sti.SyphTest(product=make_syph_dx(), test_prob_data=0.2, start=2000)
    sim = ss.Sim(n_agents=n_agents, start=2000, stop=2005, dt=ss.months(1),
                 diseases=sti.Syphilis(), networks=sti.StructuredSexual(),
                 interventions=test, demographics=[ss.Births(), ss.Deaths()], verbose=0)
    sim.run()

    assert sim.results.syphtest.new_tests.sum() > 0, 'Expected some agents to be tested with default eligibility'
    return sim


def test_syphtest_stratified_dataframe():
    """ SyphTest accepts a long-format DataFrame stratified by risk_group/sex/sw """
    sc.heading('Testing SyphTest with a stratified test_prob_data DataFrame...')

    test_prob_data = pd.DataFrame([
        {'year': 2000, 'risk_group': rg, 'sex': sex, 'sw': sw, 'symp_test_prob': 0.15}
        for rg in [0, 1, 2] for sex in ['female', 'male'] for sw in [0, 1]
    ] + [
        {'year': 2024, 'risk_group': rg, 'sex': sex, 'sw': sw, 'symp_test_prob': 0.30}
        for rg in [0, 1, 2] for sex in ['female', 'male'] for sw in [0, 1]
    ])

    test = sti.SyphTest(product=make_syph_dx(), test_prob_data=test_prob_data, start=2000)
    sim = ss.Sim(n_agents=n_agents, start=2000, stop=2010, dt=ss.months(1),
                 diseases=sti.Syphilis(), networks=sti.StructuredSexual(),
                 interventions=test, demographics=[ss.Births(), ss.Deaths()], verbose=0)
    sim.run()

    assert sim.results.syphtest.new_tests.sum() > 0, 'Expected some agents to be tested from the stratified DataFrame'
    return sim


if __name__ == '__main__':
    test_stitest_default_eligibility()
    test_syphtest_stratified_dataframe()
