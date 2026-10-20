"""
Testing and treatment of the bacterial STIs: STIDx, STITest, STITreatment, SymptomaticTesting,
SyndromicManagement, and GonorrheaTreatment with resistance
"""
import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

kw = dict(n_agents=5000, start=2000, stop=2010, verbose=0)
betas = dict(ng=dict(beta_m2f=0.4), ct=dict(beta_m2f=0.3), tv=dict(beta_m2f=0.3))


def make_sim(interventions=None, **kwargs):
    return sti.Sim(diseases=['ng', 'ct', 'tv'], **betas, init_prev=0.05, interventions=interventions, **kw, **kwargs)


def test_symptomatic_testing():
    sc.heading('Testing symptomatic testing for NG, CT, and TV, with a treatment for each')
    txs = [sti.STITreatment(diseases=d, name=f'{d}_tx') for d in ['ng', 'ct', 'tv']]
    panel = sti.SymptomaticTesting(diseases=['ng', 'ct', 'tv'], treatments=txs, name='panel')
    base, sim = ss.parallel(make_sim(), make_sim([panel] + txs)).sims
    r = sim.results
    seekers = sum(r[d].new_care_seekers.sum() for d in ['ng', 'ct', 'tv'])
    assert 0 < r.panel.delivered.sum() <= seekers # Each care seeker is tested once, whatever they sought care for
    for d in ['ng', 'ct', 'tv']:
        tx = r[f'{d}_tx']
        assert np.isclose(tx.cured.sum() + tx.failed.sum() + tx.overtreated.sum(), tx.delivered.sum())
        assert np.isclose(tx.cured.sum()/(tx.cured.sum() + tx.failed.sum()), 0.9, atol=0.03)
        assert tx.overtreated.sum() > 0 # False positives, from the specificity
        assert np.isclose(r[d].cure.sum(), tx.cured.sum())
        assert np.array_equal(tx.new_treated_success, tx.cured) # Results as STIsim 3 named them
        assert r[d].prevalence[-1] < base.results[d].prevalence[-1], f'{d}: treatment should lower prevalence'
    return sim


def test_syndromic():
    sc.heading('Testing syndromic management')
    txs = [sti.STITreatment(diseases='ng', name='ng_tx'), sti.STITreatment(diseases='ct', name='ct_tx'), sti.STITreatment(diseases='tv', name='mtz')]
    synd = sti.SyndromicManagement(diseases=['ng', 'ct', 'tv'], treatments=txs, outcome_tx_map=dict(all3=txs, ngct=txs[:2], mtnz=txs[2:], none=[]), name='synd')
    sim = make_sim([synd] + txs).run()
    r = sim.results
    n = [r.synd[f'new_tx{i}'].sum() for i in range(4)]
    assert np.isclose(sum(n), r.synd.new_care_seekers.sum())
    assert np.isclose(r.ng_tx.delivered.sum(), n[2] + n[3]) and np.isclose(r.mtz.delivered.sum(), n[1] + n[3]) # Referred by outcome
    assert r.ct_tx.overtreated.sum() > r.ct_tx.cured.sum()*0.2 # Treated for CT because of NG or TV symptoms
    assert np.isclose(r.synd.new_care_seekers_f.sum() + r.synd.new_care_seekers_m.sum(), r.synd.new_care_seekers.sum())
    return sim


def test_screening():
    sc.heading('Testing screening with STITest, to treatment, with a capacity')
    df = pd.DataFrame(dict(state=['susceptible', 'infected'], p_positive=[0.02, 0.9]))
    tx = sti.STITreatment(diseases='ct', name='ct_tx')
    young_women = ss.ppl.female & (ss.ppl.age >= 15) & (ss.ppl.age < 25)
    screen = sti.STITest(sti.STIDx(df, disease='ct'), test_prob_data=0.5, eligibility=young_women, name='ct_screen')
    sim = make_sim([screen]).run()
    r = sim.results
    assert r.ct_screen.new_tests.sum() > 0 and np.array_equal(r.ct_screen.new_tests, r.ct_screen.delivered)
    assert np.isclose(r.ct_tx.delivered.sum(), r.ct_screen.diagnosed.sum())
    assert r.ct_tx.overtreated.sum() > 0 # False positives
    yearly = r.ct_screen.delivered[1:13].sum()/r.ct_screen.eligible[1:13].mean()
    assert np.isclose(yearly, -np.log(0.5), rtol=0.1), f'Expected tests at a rate of {-np.log(0.5):.2f} per year among those eligible (a probability of 0.5 a year), not {yearly:.2f}'

    # A capacity on the treatment: at most 5 a step
    tx = sti.STITreatment(diseases='ct', max_capacity=5, name='ct_tx')
    screen = sti.STITest(sti.STIDx(df, disease='ct'), test_prob_data=0.5, eligibility=young_women, name='ct_screen')
    capped = make_sim([screen]).run()
    assert capped.results.ct_tx.delivered.max() <= 5 and capped.results.ct_tx.dropped.sum() > 0
    return sim


def test_amr():
    sc.heading('Testing gonorrhea resistance as a strain')
    ng = sti.Gonorrhea(amr=True, init_res=0.02, beta_m2f=0.4, init_prev=0.05)
    tx = sti.GonorrheaTreatment(amr=True, res_treat_eff=0.3, p_acquire=0.2, name='ng_tx')
    panel = sti.SymptomaticTesting(diseases=['ng'], treatments=[tx], name='panel')
    sim = sti.Sim(diseases=ng, interventions=[panel, tx], n_agents=5000, start=2000, stop=2015, verbose=0).run()
    r = sim.results
    res, inf = r.ng.infection.select(amr='res'), r.ng.infection
    share = lambda s: res[s].sum()/inf[s].sum()
    assert share(slice(1, 13)) < 0.1 and share(slice(-12, None)) > 0.5, 'Resistance should spread under treatment'
    for level,eff in [('sens', 0.96), ('res', 0.3)]:
        c, f = r.ng_tx.cured.select(amr=level).sum(), r.ng_tx.failed.select(amr=level).sum()
        assert np.isclose(c/(c + f), eff, atol=0.05), f'{level}: cure share {c/(c + f):.2f}, expected {eff}'
    assert r.ng_tx_fail.delivered.sum() > 0 # Resistance acquired on failure
    return sim


if __name__ == '__main__':
    sc.options(interactive=False)
    test_symptomatic_testing()
    test_syndromic()
    test_screening()
    test_amr()
