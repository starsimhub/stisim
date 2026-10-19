"""
Testing and treatment for the bacterial STIs, as Starsim products and deliveries
"""

import numpy as np
import sciris as sc
import starsim as ss
from stisim.diseases.sti import sexed

__all__ = ['STIDx', 'STITest', 'STITreatment', 'SymptomaticTesting', 'SyndromicManagement', 'care_seekers']

# The states of sti.SEIS that STIsim 3's state names cover; later rows of a test's dataframe override earlier ones
v3_states = dict(
    exposed      = ['E'],
    asymptomatic = ['E', 'P', 'A'],
    symptomatic  = ['Y', 'C'],
    seeking_care = ['C'],
    pid          = ['PID'],
    infectious   = ['P', 'A', 'Y', 'C', 'PID'],
    infected     = ['E', 'P', 'A', 'Y', 'C', 'PID'],
)


def care_seekers(diseases):
    """
    Who sought care for any of these diseases (entered C, since this point in the last step), as a function of the sim; for `eligible=`

    Args:
        diseases (list): names of `sti.SEIS` diseases
    """
    names = [getattr(d, 'name', d) for d in sc.tolist(diseases)]
    def seeking(sim):
        out = ss.uids()
        for name in names:
            out = out | sim.diseases[name].graph.states.C.enter.uids
        return out
    return seeking


def by_sex(val):
    """ A value by sex as `[women, men]` (or a scalar), as a Starsim value """
    return sexed(val) if isinstance(val, (list, tuple)) else val


def add_aliases(delivery, aliases):
    """ Copy results to the names STIsim 3 gave them """
    res = delivery.results
    for old,new in aliases.items():
        if new in res and old not in res:
            r = res[new]
            delivery.results += ss.Result(old, values=sc.dcp(r.values), scale=False, label=r.label, timevec=r.timevec, auto_plot=False, module=delivery.name)
    return


class STIDx(ss.Test):
    """
    A diagnostic test for an STI

    Either give the sensitivity and specificity (each a number, `[women, men]`, or `ss.by`,
    and the sensitivity optionally a dict by state), or, as in STIsim 3, a dataframe with
    columns `state` and `p_positive`: the row for `susceptible` gives 1 - specificity, and
    the others the sensitivity in the states they name (`exposed`, `asymptomatic`,
    `symptomatic`, `pid`, `infectious`, `infected`, or the states of `sti.SEIS`).

    Args:
        df (dataframe): the probability of a positive result by state, as in STIsim 3
        disease (str): the disease it tests for; passed on to `sti.STITest`
        sensitivity (float/list/dict): the probability of a positive result for the infected
        specificity (float/list): the probability of a negative result for the uninfected
        kwargs (dict): passed to `ss.Test` (e.g. `delay=`, `cost=`, `name=`)

    Examples:
        ```python
        naat = sti.STIDx(sensitivity=0.95, specificity=0.99, disease='ng')
        df = pd.DataFrame(dict(state=['susceptible', 'infected'], p_positive=[0.02, 0.9]))
        dx = sti.STIDx(df, disease='ct')
        ```
    """
    def __init__(self, df=None, disease=None, sensitivity=1.0, specificity=1.0, **kwargs):
        if df is not None:
            sensitivity = {}
            for state,p in zip(df.state, df.p_positive):
                if state == 'susceptible':
                    specificity = 1 - p
                else:
                    for s in v3_states.get(state, [state]):
                        sensitivity[s] = p
        if isinstance(sensitivity, dict):
            sensitivity = {k:by_sex(v) for k,v in sensitivity.items()}
        super().__init__(sensitivity=by_sex(sensitivity), specificity=by_sex(specificity), **kwargs)
        self.for_disease = disease
        return


class STITest(ss.screen):
    """
    Testing for an STI: each step, each eligible person is tested with a probability

    By default (`dt_scale=True`), `test_prob_data` is an annual probability of testing,
    applied as a rate at every step; with `dt_scale=False` it is a probability per step.
    What happens to people who test positive is `on_positive=` (e.g. an `sti.STITreatment`).

    Args:
        product (`sti.STIDx`): the test
        test_prob_data (float/array): the probability of testing; an array gives one per year in `years` (the nearest year is used)
        years (array): the years of `test_prob_data`; also the start and stop if those are not given
        start (float): the year testing starts
        stop (float): the year testing stops
        eligibility (Cond/func): who can be tested: a condition, or a function of the sim returning a boolean array or UIDs (default: everyone)
        rel_test (float): a multiplier on the probability of testing
        dt_scale (bool): whether `test_prob_data` is an annual probability (True) or a probability per step
        disease (str): the disease tested for (default: the product's)
        kwargs (dict): passed to `ss.screen` (e.g. `on_positive=`, `capacity=`, `name=`)

    Results: those of `ss.screen` (`delivered`, `diagnosed`, ...), and `new_tests` and `new_diagnoses` as in STIsim 3.

    Examples:
        ```python
        tx = sti.STITreatment(diseases='ng', name='ng_tx')
        test = sti.STITest(sti.STIDx(sensitivity=0.95, disease='ng'), test_prob_data=0.2, eligibility=ss.people.female, on_positive=tx)
        ```
    """
    def __init__(self, product=None, test_prob_data=None, years=None, start=None, stop=None, eligibility=None,
                 rel_test=1.0, dt_scale=True, disease=None, **kwargs):
        if years is not None and start is not None:
            errormsg = 'Provide either years or start, not both.'
            raise ValueError(errormsg)
        if years is not None:
            start = years[0]
            stop = sc.ifelse(stop, years[-1])
        if product is None:
            product = STIDx(disease=disease)
        self.test_prob_data = 1.0 if test_prob_data is None else test_prob_data
        self.years = years
        self.rel_test = rel_test
        self.dt_scale = dt_scale
        super().__init__(product, start=start, stop=stop, eligible=eligibility, coverage=self.test_prob, repeat=True,
                         disease=sc.ifelse(disease, product.for_disease if isinstance(product, STIDx) else None), **kwargs)
        return

    def test_prob(self, sim):
        """ The probability of testing now: per year (as a rate), or per step """
        data = self.test_prob_data
        if sc.isnumber(data):
            p = data
        else:
            p = data[sc.findnearest(self.years, sim.t.now('year'))]
        p = np.clip(p*self.rel_test, 0, 1)
        return ss.probperyear(p) if self.dt_scale else p

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_tests='delivered', new_diagnoses='diagnosed'))
        return


class STITreatment(ss.treat):
    """
    Treatment of one or more STIs: cures with probability `treat_eff`, along each disease's `cure`

    People are treated if they are referred (by a test's `on_positive=`, or by `refer()`,
    which `sti.SymptomaticTesting` and `sti.SyndromicManagement` use), or, with
    `eligibility=`, at each step that they are eligible. A treatment with no `eligibility`
    only treats people it is referred.

    Args:
        diseases (str/list): the diseases it treats
        treat_eff (float/list): the probability of cure (or `[women, men]`, or `ss.by`)
        treat_prob (float): the probability that someone referred or eligible accepts treatment
        eligibility (Cond/func): who is treated at each step; if None, only those referred
        max_capacity (float/ss.Capacity): the most that can be treated (see `ss.Delivery` `capacity=`)
        years (list): the first and last years of treatment
        on_fail (Transition/Delivery/func): what happens to people not cured (see `ss.Treatment`)
        product (`ss.Treatment`): the treatment itself, instead of one made from `diseases` and `treat_eff`
        kwargs (dict): passed to `ss.treat` (e.g. `queue=`, `priority=`, `name=`)

    Results: those of `ss.treat` (`delivered`, `cured`, `failed`, `overtreated`, by disease
    with `.select(disease=...)`), and `new_treated`, `new_treated_success`,
    `new_treated_failure`, and `new_treated_unnecessary` as in STIsim 3 (summed over diseases).

    Examples:
        ```python
        ng_tx = sti.STITreatment(diseases='ng', name='ng_tx')
        mtz = sti.STITreatment(diseases=['tv', 'bv'], treat_eff=0.95, name='metronidazole')
        ```
    """
    def __init__(self, diseases=None, treat_eff=0.9, treat_prob=1.0, eligibility=None, max_capacity=None, years=None,
                 on_fail=None, product=None, pars=None, name=None, label=None, **kwargs):
        pars = sc.mergedicts(pars)
        treat_eff = by_sex(sexed(pars.pop('treat_eff', treat_eff)))
        treat_prob = sexed(pars.pop('treat_prob', treat_prob))
        if pars:
            errormsg = f'Unknown parameters for STITreatment: {sc.strjoin(pars.keys())}'
            raise ValueError(errormsg)
        self.treats = [getattr(d, 'name', d) for d in sc.tolist(diseases)]
        if product is None:
            if not self.treats:
                errormsg = 'STITreatment needs diseases= (or a product=)'
                raise ValueError(errormsg)
            product = ss.Treatment(moves={f'{d}.cure':treat_eff for d in self.treats}, on_fail=on_fail, name=f'{name or "_".join(self.treats)}_drug')
        start, stop = (years[0], years[-1]) if years is not None else (kwargs.pop('start', None), kwargs.pop('stop', None))
        super().__init__(product, eligible=eligibility, coverage=None if treat_prob == 1 else treat_prob, capacity=max_capacity,
                         start=start, stop=stop, repeat=True, name=name, label=label, **kwargs)
        if eligibility is None and self.on is None and self.deliver is None:
            self.triggered = True # Only those referred to it
        return

    def refer(self, uids):
        """ Offer treatment now to these people (subject to treat_prob and capacity); agents only """
        self.prepare()
        self.receive(self.selector.to_selection(ss.uids(uids)))
        return

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_treated='delivered', new_treated_success='cured', new_treated_failure='failed', new_treated_unnecessary='overtreated'))
        return


def disease_name(d):
    return getattr(d, 'name', d)


class SymptomaticTesting(ss.screen):
    """
    Etiological testing of people seeking care, for several STIs at one visit

    Each person tested is tested for every disease (one `ss.screen` per disease, at the same
    visit: see `then=`), and those who test positive for a disease are referred to its
    treatment.

    Args:
        diseases (list): the diseases tested for (names or modules)
        treatments (list): `sti.STITreatment`s; the treatment of each disease is the one that treats it, unless given by `disease_treatment_map`
        disease_treatment_map (dict): the treatment of each disease, by name
        sens (dict): the sensitivity for each disease, as `[women, men]`
        spec (dict): the specificity for each disease, as `[women, men]`
        eligibility (Cond/func): who is tested (default: those who sought care for any of the diseases in this step)
        kwargs (dict): passed to `ss.screen` (e.g. `start=`, `stop=`, `capacity=`, `name=`)

    Examples:
        ```python
        ng_tx, ct_tx, tv_tx = [sti.STITreatment(diseases=d, name=f'{d}_tx') for d in ['ng', 'ct', 'tv']]
        panel = sti.SymptomaticTesting(diseases=['ng', 'ct', 'tv'], treatments=[ng_tx, ct_tx, tv_tx])
        sim = sti.Sim(diseases=['ng', 'ct', 'tv'], interventions=[panel, ng_tx, ct_tx, tv_tx])
        ```
    """
    def __init__(self, diseases=None, treatments=None, disease_treatment_map=None, sens=None, spec=None, eligibility=None,
                 negative_treatments=None, p_mtnz=None, years=None, name=None, label=None, **kwargs):
        if negative_treatments or p_mtnz:
            errormsg = 'Treating people who test negative for every disease (negative_treatments, p_mtnz) is not supported yet in Starsim 4'
            raise NotImplementedError(errormsg)
        names = [disease_name(d) for d in sc.tolist(diseases)]
        if not names:
            errormsg = 'SymptomaticTesting needs diseases='
            raise ValueError(errormsg)
        sens = sc.mergedicts({d:[0.98, 0.98] for d in names}, sens)
        spec = sc.mergedicts({d:[0.95, 0.95] for d in names}, spec)
        txmap = {d:t for t in sc.tolist(treatments) for d in getattr(t, 'treats', [])}
        txmap.update(disease_treatment_map or {})
        self.treatments = sc.tolist(treatments)
        self.disease_names = names
        if years is not None:
            kwargs.setdefault('start', years[0])
            kwargs.setdefault('stop', years[-1])
        tests = {d:STIDx(sensitivity=sens[d], specificity=spec[d], name=f'{name or "symptomatic"}_{d}_dx') for d in names}
        others = [ss.screen(tests[d], disease=d, on_positive=txmap.get(d), name=f'{name}_{d}' if name else None) for d in names[1:]]
        super().__init__(tests[names[0]], disease=names[0], eligible=sc.ifelse(eligibility, care_seekers(names)), repeat=True,
                         on_positive=txmap.get(names[0]), then=others or None, name=name, label=label, **kwargs)
        return


class SyndromicManagement(ss.Intervention):
    """
    Syndromic management of vaginal and urethral discharge (agents only)

    Each person seen is given one of four outcomes, by sex and by whether they have
    symptoms of a cervical infection (NG or CT by default):

        all3: treated for the cervical infections and for TV/BV
        ngct: treated for the cervical infections only
        mtnz: treated with metronidazole only (TV/BV)
        none: not treated

    and is referred to the treatments of their outcome (`outcome_tx_map`). The outcome
    probabilities are `tx_mix_cerv` (women with cervical symptoms, and all men) and
    `tx_mix_noncerv` (other women), each a dict of `[women, men]` by outcome.

    Args:
        treatments (list): `sti.STITreatment`s; by default, [cervical, other, ...] in outcome_tx_map
        diseases (list): the diseases whose care seekers are seen (names or modules)
        cervical_diseases (list): the diseases whose symptoms count as cervical (default: NG and CT)
        outcome_tx_map (dict): the treatments of each outcome
        tx_mix_cerv (dict): the probability of each outcome for women with cervical symptoms (and men)
        tx_mix_noncerv (dict): the probability of each outcome for other women
        eligibility (func): who is seen (default: those who sought care for any of the diseases in this step)
        start (float): the year it starts
        stop (float): the year it stops

    Results: `new_care_seekers`, and `new_tx0` to `new_tx3` (none, metronidazole only, NG/CT only, all three), with `_f` and `_m`.
    """
    def __init__(self, pars=None, treatments=None, diseases=None, cervical_diseases=None, outcome_tx_map=None,
                 eligibility=None, years=None, start=None, stop=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            tx_mix_cerv = dict(all3=[0.50, 0.05], ngct=[0.20, 0.80], mtnz=[0.20, 0.00], none=[0.10, 0.15]),
            tx_mix_noncerv = dict(all3=[0.40, 0.05], ngct=[0.10, 0.80], mtnz=[0.20, 0.00], none=[0.30, 0.15]),
            tx_cerv_f = ss.choice(a=4),
            tx_cerv_m = ss.choice(a=4),
            tx_noncerv_f = ss.choice(a=4),
        )
        self.update_pars(pars, **kwargs)
        self.treatments = sc.tolist(treatments)
        self.diseases = [disease_name(d) for d in sc.tolist(diseases)]
        self.cervical = [disease_name(d) for d in sc.tolist(cervical_diseases)] or [d for d in self.diseases if d in ['ng', 'ct']]
        if outcome_tx_map is None:
            t = self.treatments
            outcome_tx_map = dict(all3=t, ngct=t[:2], mtnz=t[1:2], none=[])
        self.outcome_tx_map = outcome_tx_map
        for tx in self.treatments: # Treated only when referred
            tx.triggered = True
        self.eligibility = sc.ifelse(eligibility, care_seekers(self.diseases))
        if years is not None:
            start, stop = sc.ifelse(start, years[0]), sc.ifelse(stop, years[-1])
        self.start, self.stop = start, stop
        self.outcomes = None
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        p = self.pars
        self.pars.tx_cerv_f.set(p=[v[0] for v in p.tx_mix_cerv.values()])
        self.pars.tx_cerv_m.set(p=[v[1] for v in p.tx_mix_cerv.values()])
        self.pars.tx_noncerv_f.set(p=[v[0] for v in p.tx_mix_noncerv.values()])
        return

    def init_results(self):
        super().init_results()
        results = []
        for sk in ['', '_f', '_m']:
            sl = f' ({sk[1:].upper()})' if sk else ''
            results += [ss.Result(f'new_care_seekers{sk}', dtype=int, label=f'Care seekers{sl}')]
            results += [ss.Result(f'new_tx{i}{sk}', dtype=int, label=f'{lbl}{sl}') for i,lbl in enumerate(['No treatment', '1 treatment', '2 treatments', '3 treatments'])]
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ppl = sim.people
        self.outcomes = None
        year = sim.t.now('year')
        if (self.start is not None and year < self.start) or (self.stop is not None and year >= self.stop):
            return
        uids = self.check_eligibility()
        if not len(uids):
            return
        f_uids = uids[ppl.female[uids]]
        m_uids = uids[~ppl.female[uids]]
        cerv = np.zeros(len(f_uids), dtype=bool)
        for d in self.cervical:
            cerv |= sim.diseases[d].symptomatic[f_uids]
        f_cerv, f_noncerv = f_uids[cerv], f_uids[~cerv]
        ofc = self.pars.tx_cerv_f.rvs(f_cerv)
        ofn = self.pars.tx_noncerv_f.rvs(f_noncerv)
        om = self.pars.tx_cerv_m.rvs(m_uids)
        self.outcomes = {key:f_cerv[ofc == i] | f_noncerv[ofn == i] | m_uids[om == i] for i,key in enumerate(['all3', 'ngct', 'mtnz', 'none'])}
        for key,txs in self.outcome_tx_map.items():
            for tx in txs:
                if len(self.outcomes[key]):
                    tx.refer(self.outcomes[key])

        ti = self.ti
        for sk,sel in [('', ppl.alive), ('_f', ppl.female), ('_m', ~ppl.female)]:
            self.results[f'new_care_seekers{sk}'][ti] += np.count_nonzero(sel[uids])
            for i,key in enumerate(['none', 'mtnz', 'ngct', 'all3']):
                self.results[f'new_tx{i}{sk}'][ti] += np.count_nonzero(sel[self.outcomes[key]])
        return
