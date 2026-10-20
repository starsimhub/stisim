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
    Who sought care for any of these diseases (entered C, since this point in the last step), as a condition; for `where=`

    Args:
        diseases (list): names of `sti.SEIS` diseases
    """
    names = [getattr(d, 'name', d) for d in sc.tolist(diseases)]
    return ss.any(ss.state(f'{name}.C').enter for name in names)


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


class STIDx(ss.Diagnostic):
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
        kwargs (dict): passed to `ss.Diagnostic` (e.g. `delay=`, `cost=`, `name=`)

    Examples:
        ```python
        naat = sti.STIDx(sensitivity=0.95, specificity=0.99, disease='ng')
        df = pd.DataFrame(dict(state=['susceptible', 'infected'], p_positive=[0.02, 0.9]))
        dx = sti.STIDx(df, disease='ct')
        ```
    """
    state_map = v3_states # The states each of STIsim 3's state names covers

    def __init__(self, df=None, disease=None, sensitivity=1.0, specificity=1.0, **kwargs):
        if df is not None:
            sensitivity = {}
            for state,p in zip(df.state, df.p_positive):
                if state == 'susceptible':
                    specificity = 1 - p
                else:
                    if state in ['trep', 'nontrep']:
                        errormsg = f'A test keyed by "{state}" is not supported: give the states it detects, or detects= (e.g. the rpr dimension of sti.Syphilis)'
                        raise NotImplementedError(errormsg)
                    for s in self.state_map.get(state, [state]):
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
    What happens to people who test positive is a step of its own, triggered by the
    result: e.g. `sti.STITreatment(diseases='ng', who=test.positive)`.

    Args:
        product (`sti.STIDx`): the test
        test_prob_data (float/array): the probability of testing; an array gives one per year in `years` (the nearest year is used)
        years (array): the years of `test_prob_data`; also the start and stop if those are not given
        start (float): the year testing starts
        stop (float): the year testing stops
        eligibility (Cond/func): who can be tested: a condition, or a function of the sim (default: everyone)
        rel_test (float): a multiplier on the probability of testing
        dt_scale (bool): whether `test_prob_data` is an annual probability (True) or a probability per step
        disease (str): the disease tested for (default: the product's)
        kwargs (dict): passed to `ss.screen` (e.g. `who=`, `capacity=`, `name=`)

    Results: those of `ss.screen` (`delivered`, `positive`, `negative`, ...), and `new_tests` and `new_diagnoses` as in STIsim 3.

    Examples:
        ```python
        test = sti.STITest(sti.STIDx(sensitivity=0.95, disease='ng'), test_prob_data=0.2, eligibility=ss.ppl.female, name='ng_test')
        tx = sti.STITreatment(diseases='ng', who=test.positive, name='ng_tx')
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
        super().__init__(product, start=start, stop=stop, where=eligibility, coverage=self.test_prob,
                         disease=sc.ifelse(disease, product.for_disease if isinstance(product, STIDx) else None), **kwargs)
        return

    def test_prob(self, sim):
        """ The probability of testing now, as a typed probability (per year, or per step), so it is offered at every step """
        data = self.test_prob_data
        if sc.isnumber(data):
            p = data
        else:
            p = data[sc.findnearest(self.years, sim.t.now('year'))]
        p = np.clip(p*self.rel_test, 0, 1)
        return ss.probperyear(p) if self.dt_scale else ss.prob(p, sim.t.dt)

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_tests='delivered', new_diagnoses='positive'))
        return


class STITreatment(ss.treat):
    """
    Treatment of one or more STIs: cures with probability `treat_eff`, along each disease's `cure`

    Who is treated is `who=`: usually the outcome of a test (`who=test.positive`) or of
    syndromic management (`who=synd.all3 | synd.ngct`), or `eligibility=`, a condition
    (treated once in each spell of meeting it).

    Args:
        diseases (str/list): the diseases it treats
        treat_eff (float/list): the probability of cure (or `[women, men]`, or `ss.by`)
        treat_prob (float): the probability that someone referred or eligible accepts treatment
        eligibility (Cond/func): who is treated, as a condition (`where=` of the delivery)
        max_capacity (float): the most that can be treated (`capacity=` of the delivery)
        years (list): the first and last years of treatment
        product (`ss.Treatment`): the treatment itself, instead of one made from `diseases` and `treat_eff`
        kwargs (dict): passed to `ss.treat` (e.g. `who=`, `queue=`, `priority=`, `uses=`, `name=`)

    Results: those of `ss.treat` (`delivered`, `cured`, `failed`, `unneeded`, by disease
    with `.select(disease=...)`), and `new_treated`, `new_treated_success`,
    `new_treated_failure`, and `new_treated_unnecessary` as in STIsim 3 (summed over diseases).

    What happens to treatment failures is a step of its own (`who=tx.failed`), or a
    transition of the disease triggered by the outcome (see `sti.GonorrheaTreatment`).

    Examples:
        ```python
        ng_tx = sti.STITreatment(diseases='ng', who=test.positive, name='ng_tx')
        mtz = sti.STITreatment(diseases=['tv', 'bv'], treat_eff=0.95, who=synd.mtnz | synd.all3, name='metronidazole')
        ```
    """
    def __init__(self, diseases=None, treat_eff=0.9, treat_prob=1.0, eligibility=None, max_capacity=None, years=None,
                 product=None, pars=None, name=None, label=None, **kwargs):
        pars = sc.mergedicts(pars)
        treat_eff = by_sex(sexed(pars.pop('treat_eff', treat_eff)))
        treat_prob = sexed(pars.pop('treat_prob', treat_prob))
        if pars:
            errormsg = f'Unknown parameters for STITreatment: {sc.strjoin(pars.keys())}'
            raise ValueError(errormsg)
        if eligibility is None and kwargs.get('who') is None and kwargs.get('on') is None and kwargs.get('deliver') is None:
            errormsg = 'STITreatment needs to know who is treated: who= (e.g. who=test.positive) or eligibility= (a condition)'
            raise ValueError(errormsg)
        self.treats = [getattr(d, 'name', d) for d in sc.tolist(diseases)]
        if product is None:
            if not self.treats:
                errormsg = 'STITreatment needs diseases= (or a product=)'
                raise ValueError(errormsg)
            product = ss.Treatment(moves={f'{d}.cure':treat_eff for d in self.treats}, name=f'{name or "_".join(self.treats)}_drug')
        start, stop = (years[0], years[-1]) if years is not None else (kwargs.pop('start', None), kwargs.pop('stop', None))
        super().__init__(product, where=eligibility, coverage=None if treat_prob == 1 else treat_prob, capacity=max_capacity,
                         start=start, stop=stop, name=name, label=label, **kwargs)
        return

    def finalize_results(self):
        super().finalize_results()
        add_aliases(self, dict(new_treated='delivered', new_treated_success='cured', new_treated_failure='failed', new_treated_unnecessary='unneeded'))
        return


def disease_name(d):
    return getattr(d, 'name', d)


class SymptomaticTesting(ss.cascade):
    """
    Etiological testing of people seeking care, for several STIs at one visit

    A cascade of a visit (`<name>_visit`, by default for those who sought care for any of
    the diseases on the step) and a test for each disease at it (`<name>_<disease>`). The
    treatment of those who test positive for a disease is a step of its own, triggered by
    `positive(disease)`.

    Args:
        diseases (list): the diseases tested for (names or modules)
        sens (dict): the sensitivity for each disease, as `[women, men]`
        spec (dict): the specificity for each disease, as `[women, men]`
        eligibility (Cond/event): who attends the visit (default: those who sought care for any of the diseases)
        years (list): the first and last years of testing
        kwargs (dict): passed to the visit (`ss.visit`, e.g. `start=`, `stop=`, `coverage=`, `capacity=`)

    Examples:
        ```python
        panel = sti.SymptomaticTesting(diseases=['ng', 'ct', 'tv'], name='panel')
        txs = [sti.STITreatment(diseases=d, who=panel.positive(d), name=f'{d}_tx') for d in ['ng', 'ct', 'tv']]
        sim = sti.Sim(diseases=['ng', 'ct', 'tv'], interventions=[panel, *txs])
        ```
    """
    def __init__(self, diseases=None, sens=None, spec=None, eligibility=None, negative_treatments=None, p_mtnz=None,
                 years=None, name='symptomatic_testing', label=None, **kwargs):
        if negative_treatments or p_mtnz:
            errormsg = 'Treating people who test negative for every disease (negative_treatments, p_mtnz) is not supported yet in Starsim 4'
            raise NotImplementedError(errormsg)
        if 'treatments' in kwargs or 'disease_treatment_map' in kwargs:
            errormsg = 'The treatments of SymptomaticTesting are steps of their own, triggered by its results: sti.STITreatment(diseases=d, who=panel.positive(d))'
            raise TypeError(errormsg)
        names = [disease_name(d) for d in sc.tolist(diseases)]
        if not names:
            errormsg = 'SymptomaticTesting needs diseases='
            raise ValueError(errormsg)
        sens = sc.mergedicts({d:[0.98, 0.98] for d in names}, sens)
        spec = sc.mergedicts({d:[0.95, 0.95] for d in names}, spec)
        if years is not None:
            kwargs.setdefault('start', years[0])
            kwargs.setdefault('stop', years[-1])
        visit = ss.visit(who=sc.ifelse(eligibility, care_seekers(names)), name=f'{name}_visit', **kwargs)
        screens = {d:ss.screen(STIDx(sensitivity=sens[d], specificity=spec[d], name=f'{name}_{d}_dx'), who=visit, disease=d, name=f'{name}_{d}') for d in names}
        super().__init__(visit, *screens.values(), name=name, label=label)
        self.disease_names = names
        self.visit = visit
        self.screens = screens
        return

    def positive(self, disease):
        """ The people who tested positive for this disease, as an event that triggers its treatment """
        return self.screens[disease_name(disease)].positive


class SyndromicManagement(ss.Intervention):
    """
    Syndromic management of vaginal and urethral discharge (agents only)

    Each person seen is given one of four outcomes, by sex and by whether they have
    symptoms of a cervical infection (NG or CT by default), each an event that the
    treatments are triggered by:

        all3: treated for the cervical infections and for TV/BV
        ngct: treated for the cervical infections only
        mtnz: treated with metronidazole only (TV/BV)
        none: not treated

    The outcome probabilities are `tx_mix_cerv` (women with cervical symptoms, and all
    men) and `tx_mix_noncerv` (other women), each a dict of `[women, men]` by outcome.
    List it before the treatments it triggers.

    Args:
        diseases (list): the diseases whose care seekers are seen (names or modules)
        cervical_diseases (list): the diseases whose symptoms count as cervical (default: NG and CT)
        tx_mix_cerv (dict): the probability of each outcome for women with cervical symptoms (and men)
        tx_mix_noncerv (dict): the probability of each outcome for other women
        eligibility (func): who is seen, as a function of the sim returning UIDs (default: those who sought care for any of the diseases on the step)
        start (float): the year it starts
        stop (float): the year it stops

    Results: `new_care_seekers`, and `new_tx0` to `new_tx3` (none, metronidazole only, NG/CT only, all three), with `_f` and `_m`.

    Examples:
        ```python
        synd = sti.SyndromicManagement(diseases=['ng', 'ct', 'tv'], name='synd')
        ng_tx = sti.STITreatment(diseases='ng', who=synd.all3 | synd.ngct, name='ng_tx')
        ct_tx = sti.STITreatment(diseases='ct', who=synd.all3 | synd.ngct, name='ct_tx')
        mtz = sti.STITreatment(diseases='tv', who=synd.all3 | synd.mtnz, name='mtz')
        sim = sti.Sim(diseases=['ng', 'ct', 'tv'], interventions=[synd, ng_tx, ct_tx, mtz])
        ```
    """
    outcome_keys = ['all3', 'ngct', 'mtnz', 'none']

    def __init__(self, pars=None, diseases=None, cervical_diseases=None, eligibility=None, years=None, start=None, stop=None,
                 name=None, label=None, **kwargs):
        if 'treatments' in kwargs or 'outcome_tx_map' in kwargs:
            errormsg = 'The treatments of SyndromicManagement are steps of their own, triggered by its outcomes: e.g. sti.STITreatment(diseases="ng", who=synd.all3 | synd.ngct)'
            raise TypeError(errormsg)
        super().__init__(name=name, label=label)
        self.define_pars(
            tx_mix_cerv = dict(all3=[0.50, 0.05], ngct=[0.20, 0.80], mtnz=[0.20, 0.00], none=[0.10, 0.15]),
            tx_mix_noncerv = dict(all3=[0.40, 0.05], ngct=[0.10, 0.80], mtnz=[0.20, 0.00], none=[0.30, 0.15]),
            tx_cerv_f = ss.choice(a=4),
            tx_cerv_m = ss.choice(a=4),
            tx_noncerv_f = ss.choice(a=4),
        )
        self.update_pars(pars, **kwargs)
        self.define_states(*[ss.Event(key, label=f'Syndromic management outcome: {key}') for key in self.outcome_keys])
        self.diseases = [disease_name(d) for d in sc.tolist(diseases)]
        self.cervical = [disease_name(d) for d in sc.tolist(cervical_diseases)] or [d for d in self.diseases if d in ['ng', 'ct']]
        self.eligibility = eligibility
        if years is not None:
            start, stop = sc.ifelse(start, years[0]), sc.ifelse(stop, years[-1])
        self.start, self.stop = start, stop
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        p = self.pars
        self.pars.tx_cerv_f.set(p=[p.tx_mix_cerv[k][0] for k in self.outcome_keys])
        self.pars.tx_cerv_m.set(p=[p.tx_mix_cerv[k][1] for k in self.outcome_keys])
        self.pars.tx_noncerv_f.set(p=[p.tx_mix_noncerv[k][0] for k in self.outcome_keys])
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

    def seen(self):
        """ Who is seen on this step: those who sought care for any of the diseases """
        if self.eligibility is not None:
            return ss.uids(self.eligibility(self.sim))
        return ss.uids.concatenate([self.sim.diseases[name].graph.states.C.enter.uids for name in self.diseases]).unique()

    def step(self):
        sim = self.sim
        ppl = sim.people
        year = sim.t.now('year')
        if (self.start is not None and year < self.start) or (self.stop is not None and year >= self.stop):
            return
        uids = self.seen()
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
        outcomes = {key:f_cerv[ofc == i] | f_noncerv[ofn == i] | m_uids[om == i] for i,key in enumerate(self.outcome_keys)}
        for key,out in outcomes.items():
            if len(out):
                getattr(self, key).fire(out)

        ti = self.ti
        for sk,sel in [('', ppl.alive), ('_f', ppl.female), ('_m', ~ppl.female)]:
            self.results[f'new_care_seekers{sk}'][ti + 1] += np.count_nonzero(sel[uids])
            for i,key in enumerate(['none', 'mtnz', 'ngct', 'all3']):
                self.results[f'new_tx{i}{sk}'][ti + 1] += np.count_nonzero(sel[outcomes[key]])
        return
