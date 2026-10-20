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
    def interventions():
        panel = sti.SymptomaticTesting(diseases=['ng', 'ct', 'tv'], name='panel')
        txs = [sti.STITreatment(diseases=d, who=panel.positive(d), name=f'{d}_tx') for d in ['ng', 'ct', 'tv']]
        return [panel, *txs]
    base, sim = ss.parallel(make_sim(), make_sim(interventions())).sims
    r = sim.results
    seekers = sum(r[d].new_care_seekers.sum() for d in ['ng', 'ct', 'tv'])
    assert 0 < r.panel_visit.delivered.sum() <= seekers # Each care seeker is tested once, whatever they sought care for
    for d in ['ng', 'ct', 'tv']:
        assert np.isclose(r[f'panel_{d}'].delivered.sum(), r.panel_visit.delivered.sum()) # Tested for each at the visit
        tx = r[f'{d}_tx']
        assert np.isclose(tx.delivered.sum(), r[f'panel_{d}'].positive.sum())
        assert np.isclose(tx.cured.sum() + tx.failed.sum() + tx.unneeded.sum(), tx.delivered.sum())
        assert np.isclose(tx.cured.sum()/(tx.cured.sum() + tx.failed.sum()), 0.9, atol=0.03)
        assert tx.unneeded.sum() > 0 # False positives, from the specificity
        assert np.isclose(r[d].cure.sum(), tx.cured.sum())
        assert np.array_equal(tx.new_treated_success, tx.cured) # Results as STIsim 3 named them
        assert r[d].prevalence[-1] < base.results[d].prevalence[-1], f'{d}: treatment should lower prevalence'
    return sim


def test_syndromic():
    sc.heading('Testing syndromic management')
    synd = sti.SyndromicManagement(diseases=['ng', 'ct', 'tv'], name='synd')
    txs = [sti.STITreatment(diseases='ng', who=synd.all3 | synd.ngct, name='ng_tx'),
           sti.STITreatment(diseases='ct', who=synd.all3 | synd.ngct, name='ct_tx'),
           sti.STITreatment(diseases='tv', who=synd.all3 | synd.mtnz, name='mtz')]
    sim = make_sim([synd] + txs).run()
    r = sim.results
    n = [r.synd[f'new_tx{i}'].sum() for i in range(4)]
    assert np.isclose(sum(n), r.synd.new_care_seekers.sum())
    assert np.isclose(r.ng_tx.delivered.sum(), n[2] + n[3]) and np.isclose(r.mtz.delivered.sum(), n[1] + n[3]) # Triggered by outcome
    assert r.ct_tx.unneeded.sum() > r.ct_tx.cured.sum()*0.2 # Treated for CT because of NG or TV symptoms
    assert np.isclose(r.synd.new_care_seekers_f.sum() + r.synd.new_care_seekers_m.sum(), r.synd.new_care_seekers.sum())
    return sim


def test_screening():
    sc.heading('Testing screening with STITest, to treatment, with a capacity')
    df = pd.DataFrame(dict(state=['susceptible', 'infected'], p_positive=[0.02, 0.9]))
    young_women = ss.ppl.female & (ss.ppl.age >= 15) & (ss.ppl.age < 25)
    screen = sti.STITest(sti.STIDx(df, disease='ct'), test_prob_data=0.5, eligibility=young_women, name='ct_screen')
    tx = sti.STITreatment(diseases='ct', who=screen.positive, name='ct_tx')
    sim = make_sim([screen, tx]).run()
    r = sim.results
    assert r.ct_screen.new_tests.sum() > 0 and np.array_equal(r.ct_screen.new_tests, r.ct_screen.delivered)
    assert np.isclose(r.ct_tx.delivered.sum(), r.ct_screen.positive.sum())
    assert r.ct_tx.unneeded.sum() > 0 # False positives
    yearly = r.ct_screen.delivered[1:13].sum()/r.ct_screen.eligible[1:13].mean()
    assert np.isclose(yearly, -np.log(0.5), rtol=0.1), f'Expected tests at a rate of {-np.log(0.5):.2f} per year among those eligible (a probability of 0.5 a year), not {yearly:.2f}'

    # A capacity on the treatment: at most 5 a step
    screen = sti.STITest(sti.STIDx(df, disease='ct'), test_prob_data=0.5, eligibility=young_women, name='ct_screen')
    tx = sti.STITreatment(diseases='ct', who=screen.positive, max_capacity=5, name='ct_tx')
    capped = make_sim([screen, tx]).run()
    assert capped.results.ct_tx.delivered.max() <= 5 and capped.results.ct_tx.lost_capacity.sum() > 0
    return sim


def test_amr():
    sc.heading('Testing gonorrhea resistance as a strain')
    ng = sti.Gonorrhea(amr=True, init_res=0.02, beta_m2f=0.4, init_prev=0.05)
    panel = sti.SymptomaticTesting(diseases=['ng'], name='panel')
    tx = sti.GonorrheaTreatment(amr=True, res_treat_eff=0.3, p_acquire=0.2, who=panel.positive('ng'), name='ng_tx')
    sim = sti.Sim(diseases=ng, interventions=[panel, tx, tx.acquire], n_agents=5000, start=2000, stop=2015, verbose=0).run()
    r = sim.results
    res, inf = r.ng.infection.select(amr='res'), r.ng.infection
    share = lambda s: res[s].sum()/inf[s].sum()
    assert share(slice(1, 13)) < 0.1 and share(slice(-12, None)) > 0.5, 'Resistance should spread under treatment'
    for level,eff in [('sens', 0.96), ('res', 0.3)]:
        c, f = r.ng_tx.cured.select(amr=level).sum(), r.ng_tx.failed.select(amr=level).sum()
        assert np.isclose(c/(c + f), eff, atol=0.05), f'{level}: cure share {c/(c + f):.2f}, expected {eff}'
    n_fail = r.ng_tx.failed.select(amr='sens').sum() # Resistance acquired on failure, with p_acquire: binomial, so within 3 standard deviations (a few dozen failures)
    assert abs(r.ng_tx_fail.delivered.sum() - 0.2*n_fail) <= 3*np.sqrt(n_fail*0.2*0.8) + 1, f'{r.ng_tx_fail.delivered.sum()} acquired resistance of {n_fail} failures, expected 0.2'
    return sim


def test_syphilis_care():
    sc.heading('Testing syphilis screening, treatment of mothers and fetuses, and newborn testing')
    dual = sti.SyphDx(sensitivity=dict(E=0.2, Pv=0.9, Ph=0.9, S2v=0.95, S2h=0.95, EL=0.95, LL=0.95, T=0.95), specificity=0.99, name='dual')
    test = sti.SyphTest(product=dual, test_prob_data=0.2, name='syph_test')
    anc = sti.ANCSyphTest(product=dual, test_prob_data=0.9, name='anc')
    tx = sti.SyphTx(who=test.positive | anc.positive, name='syph_tx')
    nb = sti.NewbornSyphTest(anc=anc, test_prob_data=1.0, name='nb_test')
    nb_tx = sti.NewbornTreatment(who=nb.positive, name='nb_tx')
    syph = sti.Syphilis(beta_m2f=0.3, beta_m2c=0.5, init_prev=0.1, init_latent_prev=0.1)
    kw = dict(n_agents=5000, start=2000, stop=2010, verbose=0, demographics=[ss.Pregnancy(fertility_rate=ss.freqperyear(80)), ss.Deaths(death_rate=10)])
    base, sim = ss.parallel(sti.Sim(diseases=syph, **kw), sti.Sim(diseases=syph, interventions=[test, anc, tx, tx.fetal, nb, nb_tx], **kw)).sims
    r = sim.results
    assert r.anc.delivered.sum() > 0 and r.anc.delivered.sum() <= r.pregnancy.pregnancies.sum()
    assert np.isclose(r.syph_tx.delivered.sum(), r.syph_test.positive.sum() + r.anc.positive.sum(), rtol=0.05) # Triggered by either (the same step counts once)
    assert np.isclose(r.syph.cure_early.sum() + r.syph.cure_late.sum(), r.syph_tx.cured.sum())
    assert r.syph_tx_fetus.delivered.sum() > 0 # Fetuses of cured mothers, over the maternal network
    assert np.isclose(r.syph.fetal_cure.sum(), r.syph_tx_fetus.cured.sum())
    assert r.syph.prevalence[-1] < base.results.syph.prevalence[-1]
    assert 0 < r.nb_test.delivered.sum() <= r.anc.positive.sum() and np.isclose(r.nb_tx.delivered.sum(), r.nb_test.positive.sum())
    return sim


if __name__ == '__main__':
    sc.options(interactive=False)
    test_symptomatic_testing()
    test_syndromic()
    test_screening()
    test_amr()
