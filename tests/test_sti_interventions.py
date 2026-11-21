"""
Test STI testing interventions: STITest/SyphTest eligibility and coverage input formats

Tests cover:
    - STITest.get_testers works with the default eligibility (None -> everyone), not just
      an explicit eligibility function (regression test for GH #537)
    - SyphTest accepts a DataFrame of test probabilities stratified by
      risk group/sex/sex-work status, as documented in the interventions user guide
    - SymptomaticTesting and SyndromicManagement run with their defaults and stop at stop=,
      UpdateDrugs resets treatment response, and NewbornSyphTest is not dt-scaled
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
    sim = ss.Sim(n_agents=n_agents, start=2000, stop=2003, dt=ss.months(1),
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
    sim = ss.Sim(n_agents=n_agents, start=2000, stop=2003, dt=ss.months(1),
                 diseases=sti.Syphilis(), networks=sti.StructuredSexual(),
                 interventions=test, demographics=[ss.Births(), ss.Deaths()], verbose=0)
    sim.run()

    assert sim.results.syphtest.new_tests.sum() > 0, 'Expected some agents to be tested from the stratified DataFrame'
    return sim


def test_syndromic():
    """ SymptomaticTesting, SyndromicManagement, and UpdateDrugs run with their defaults """
    sc.heading('Testing syndromic management and drug updates...')

    def symptomatic(sim):
        return (sim.diseases.ng.symptomatic | sim.diseases.ct.symptomatic | sim.diseases.tv.symptomatic).uids

    def make_sim(*intvs):
        diseases = [sti.Gonorrhea(init_prev=0.1), sti.Chlamydia(init_prev=0.1), sti.Trichomoniasis(init_prev=0.1)]
        sim = ss.Sim(n_agents=n_agents, start=2000, stop=2003, dt=ss.months(1), diseases=diseases,
                     networks=sti.StructuredSexual(), interventions=list(intvs), verbose=0)
        return sim.run()

    # SymptomaticTesting with default negative_treatments and disease_treatment_map
    txs = [sti.GonorrheaTreatment(), sti.STITreatment(diseases='ct', name='ct_tx'), sti.STITreatment(diseases='tv', name='metro')]
    sim = make_sim(sti.SymptomaticTesting(diseases=['ng', 'ct', 'tv'], eligibility=symptomatic, treatments=txs), *txs)
    assert sim.results.symptomatictesting.new_care_seekers.sum() > 0

    # SyndromicManagement stops at stop=
    txs = [sti.GonorrheaTreatment(), sti.STITreatment(diseases='tv', name='metro')]
    sim = make_sim(sti.SyndromicManagement(diseases=['ng', 'ct', 'tv'], eligibility=symptomatic, treatments=txs, stop=2002), *txs)
    seekers = sim.results.syndromicmanagement.new_care_seekers
    assert seekers.sum() > 0 and seekers[sim.t.yearvec >= 2002].sum() == 0

    # UpdateDrugs resets the treatment response once it falls below the threshold
    ng_tx = sti.GonorrheaTreatment(eligibility=lambda sim: sim.people.alive.uids, rel_treat_unneed=0.5)
    sim = make_sim(ng_tx, sti.UpdateDrugs(threshold_amr=0.3))
    assert sim.interventions.updatedrugs.change_time is not None
    return sim


def test_newborn_syph_test():
    """ NewbornSyphTest tests each scheduled newborn with probability test_prob_data, not scaled by dt """
    sc.heading('Testing NewbornSyphTest...')
    test = sti.NewbornSyphTest(test_prob_data=1.0, start=2000)
    sim = ss.Sim(n_agents=n_agents, start=2000, stop=2002, dt=ss.months(1), diseases=sti.Syphilis(),
                 networks=sti.StructuredSexual(), interventions=test, verbose=0)
    sim.init()
    test = sim.interventions.newbornsyphtest
    uids = sim.people.alive.uids
    test.schedule(uids, test.ti)
    assert len(test.get_testers(sim)) == len(uids)
    return sim


if __name__ == '__main__':
    test_stitest_default_eligibility()
    test_syphtest_stratified_dataframe()
    test_syndromic()
    test_newborn_syph_test()
