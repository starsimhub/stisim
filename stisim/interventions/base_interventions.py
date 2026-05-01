"""
Define interventions for STIsim
"""

import starsim as ss
import numpy as np
import sciris as sc
import pandas as pd
from scipy.interpolate import interp1d
from stisim.utils import count


# %% Base classes
__all__ = ["STIDx", "STITest", "SymptomaticTesting", "SyndromicMgmt", "ANCTest", "STITreatment", "PartnerNotification", "ProductMix"]


class STIDx(ss.Product):
    """
    Generic class for diagnostics with a positive/negative outcome
    Uses bernoulli sampling, so can only be used for tests with binary outcomes
    Results vary depending on agents' true underlying health state
    """
    def __init__(self, df, disease, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.disease = disease
        self.df = df
        self.health_states = df.state.unique()
        self.result_list = ['positive', 'negative']
        self.probs = {state: df.loc[df.state == state].p_positive.values[0] for state in self.health_states}

        # Store the probability of returning each possible result
        self.result_dist = ss.bernoulli(p=0)
        return

    def administer(self, sim, uids):
        """
        Administer a testing product.
        """
        p_pos = np.full(len(uids), np.nan, dtype=ss.dtypes.float)  # Default to negative
        for state in self.health_states:
            this_state = getattr(sim.diseases[self.disease], state)
            in_this_state = this_state[uids]  # Find people for which this state is true
            if in_this_state.any():
                p_pos[in_this_state] = self.probs[state]

        # Apply the test
        self.result_dist.set(p_pos)
        pos, neg = self.result_dist.split(uids)
        outcomes = {r: ss.uids() for r in self.result_list}
        outcomes['positive'] = pos
        outcomes['negative'] = neg

        return outcomes


class STITest(ss.Intervention):
    """
    Base class for STI testing.

    Controls who gets tested, how often, and with what diagnostic product.
    Each timestep, eligible agents are tested with a probability derived from
    test_prob_data. By default (dt_scale=True), test_prob_data is interpreted
    as an **annual testing probability** and converted to a per-timestep
    probability via ``ss.probperyear``. This correctly recovers the annual
    probability over a full year of timesteps (e.g. test_prob_data=0.24 with
    monthly dt gives ~2.26% per month, recovering exactly 24% per year).

    If you need to specify a per-timestep probability directly (not an annual
    probability), set dt_scale=False:
        - test_prob_data=0.5, dt_scale=False → 50% chance per timestep

    Args:
        product:        diagnostic product (e.g. STIDx, HIVDx) that determines
                        test outcomes (sensitivity, specificity). If None,
                        subclasses supply a default (e.g. HIVTest uses a perfect
                        HIVDx). For STITest, a product is required.
        test_prob_data: annual testing probability (if dt_scale=True, the default)
                        or per-timestep probability (if dt_scale=False). Accepts a
                        scalar (constant probability), an array (one value per
                        entry in ``years``), or a DataFrame (subclass-specific,
                        e.g. SymptomaticTesting accepts risk-group-stratified
                        DataFrames). Default: 1.0 (test all eligible agents per
                        year).
        years (array):  calendar years corresponding to entries in test_prob_data
                        when test_prob_data is an array. Mutually exclusive with start.
        start (float):  calendar year when the intervention activates (inclusive).
                        Defaults to the first simulation year.
        stop (float):   calendar year when the intervention deactivates (inclusive).
                        Defaults to the last simulation year.
        eligibility:    function ``f(sim) -> BoolArr`` or ``f(sim) -> UIDs`` defining
                        who can be tested. Can filter on any agent property, e.g.
                        ``lambda sim: sim.people.female & (sim.people.age > 15)``.
                        Default: all agents (subclasses may override, e.g. HIVTest
                        defaults to undiagnosed agents only).
        rel_test:       relative testing probability multiplier (default 1.0)
        dt_scale:       if True (default), interpret test_prob_data as an annual
                        probability. Set to False to interpret test_prob_data as a
                        per-timestep probability.

    States set on agents:
        tested (bool):       ever tested
        diagnosed (bool):    ever diagnosed positive
        ti_tested (float):   timestep of last test
        ti_scheduled (float): timestep of next scheduled test
        ti_positive (float): timestep of last positive result
        ti_negative (float): timestep of last negative result
        tests (float):       cumulative number of tests received
    """

    def __init__(self, pars=None, test_prob_data=None, years=None, start=None, stop=None, eligibility=None, product=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            rel_test=1,
            dt_scale=True,
        )
        self.update_pars(pars, **kwargs)

        # Years
        if years is not None and start is not None:
            errormsg = 'Provide either years or start_year, not both.'
            raise ValueError(errormsg)
        self.years = years
        self.start = start
        self.stop = stop
        if self.start is None:
            if self.years is not None:
                self.start = self.years[0]
        if self.stop is None:
            if self.years is not None:
                self.stop = self.years[-1]

        # Testing eligibility and uptake
        self.eligibility = eligibility
        if test_prob_data is None: test_prob_data = 1  # Default test everyone
        self.test_prob_data = test_prob_data
        self.test_prob = ss.bernoulli(self.make_test_prob_fn)

        # States
        self.define_states(
            ss.BoolState('tested'),
            ss.BoolState('diagnosed'),
            ss.FloatArr('ti_tested'),
            ss.FloatArr('ti_scheduled'),
            ss.FloatArr('ti_positive'),
            ss.FloatArr('ti_negative'),
            ss.FloatArr('tests', default=0),
        )

        # Products
        self.product = product
        self.outcomes = {}
        if self.product is not None:
            self.outcomes = {outcome: ss.FloatArr(f'prod_ti_{outcome}') for outcome in self.product.result_list}
            self.last_outcomes = {outcome: ss.uids() for outcome in self.product.result_list}

        return

    @property
    def state_list(self):
        """ Include products in the state list """
        out = super().state_list + list(self.outcomes.values())
        return out

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start is None: self.start = sim.t.yearvec[0]
        if self.stop is None: self.stop = sim.t.yearvec[-1]
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_diagnoses', dtype=int, label="New diagnoses", auto_plot=False),
            ss.Result('new_tests', dtype=int, label="New tests", auto_plot=False),
        )
        return

    @staticmethod
    def make_test_prob_fn(self, sim, uids):
        """ Testing probabilites over time """

        if sc.isnumber(self.test_prob_data):
            test_prob = self.test_prob_data

        elif sc.checktype(self.test_prob_data, 'arraylike'):
            year_ind = sc.findnearest(self.years, self.t.now('year'))
            test_prob = self.test_prob_data[year_ind]
        else:
            errormsg = 'Format of test_prob_data must be float or array.'
            raise ValueError(errormsg)

        # Scale and validate
        test_prob *= self.pars.rel_test
        if self.pars.dt_scale:
            test_prob = ss.probperyear(test_prob).to_prob(sim.dt)
        test_prob = np.clip(test_prob, a_min=0, a_max=1)

        return test_prob

    def get_testers(self, sim):
        """
        Find who tests by applying eligibility and coverage/uptake
        """
        # Select UIDs for testing based on eligibility and test_prob
        accept_uids = ss.uids()
        if callable(self.eligibility):
            eligible_uids = self.check_eligibility()  # Apply eligiblity - uses base class from ss.Intervention
        else:
            eligible_uids = self.eligibility
        if len(eligible_uids):
            accept_uids = self.test_prob.filter(eligible_uids)
        scheduled_uids = (self.ti_scheduled == self.ti).uids  # Add on scheduled tests
        return accept_uids | scheduled_uids

    def schedule(self, uids, ti):
        # Schedule a future test for specified UIDs at specified time indices
        self.ti_scheduled[uids] = ti

    def step(self, uids=None):
        """ Apply the testing intervention """
        sim = self.sim
        outcomes = {outcome: ss.uids() for outcome in self.product.result_list}
        self.last_outcomes = outcomes

        # Apply if within the start years
        after_start = (self.t.now('year') >= self.start)
        before_stop = (self.t.now('year') < self.stop)
        if (after_start & before_stop):

            if uids is None:
                uids = self.get_testers(sim)
                self.ti_tested[uids] = self.ti

            if len(uids):
                outcomes = self.product.administer(sim, uids)
                self.last_outcomes = outcomes
                for k in self.product.result_list:
                    self.outcomes[k][self.last_outcomes[k]] = self.ti
            else:
                return outcomes

            # If it's a binary test, set the time of positive/negative outcome
            if 'positive' in outcomes.keys():
                self.ti_negative[outcomes['negative']] = self.ti
                self.ti_positive[outcomes['positive']] = self.ti

        return outcomes

    def update_results(self):
        # Store results
        ti = self.ti
        self.results['new_diagnoses'][ti] += count(self.ti_positive == ti)
        self.results['new_tests'][ti] += count(self.ti_tested == ti)


class SymptomaticTesting(STITest):
    """
    Base class for symptomatic testing with multiple possible etiologies
    Unlike other test classes, this doesn't return positive/negative outcomes, since
    syndromic management doesn't involve reaching a diagnosis.
    Rather, the testing intervention itself contains a linked treatment intervention.
    """

    def __init__(self, pars=None, treatments=None, diseases=None, disease_treatment_map=None, negative_treatments=None,
                treat_prob_data=None, years=None, start=None, stop=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(years=years, start=start, stop=stop, eligibility=eligibility, name=name, label=label)
        self.define_pars(
            sens=dict(  # VDS: treat-all approach. UDS: treat most for NG+CT, rarely treat for TV
                ng=[0.98, 0.98],
                ct=[0.98, 0.98],
                tv=[0.98, 0.98],
            ),
            spec=dict(
                ng=[0.95, 0.95],
                ct=[0.95, 0.95],
                tv=[0.95, 0.95],
            ),
            sens_dist=ss.bernoulli(p=0),
            spec_dist=ss.bernoulli(p=0),
            p_mtnz=ss.bernoulli(p=0),
            dt_scale=False,
        )
        self.update_pars(pars, **kwargs)

        # Store treatments and diseases
        self.treatments = sc.tolist(treatments)
        self.diseases = diseases
        if disease_treatment_map is None:
            disease_treatment_map = {t.disease: t for t in self.treatments}
        self.disease_treatment_map = disease_treatment_map
        self.negative_treatments = negative_treatments

        self.define_states(
            ss.FloatArr('ti_referred'),
            ss.FloatArr('ti_dismissed'),
        )
        self.treat_prob_data = treat_prob_data
        self.treat_prob = None
        self.treated_by_uid = None
        self.n_treatments = 0

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.treat_prob_data is not None:
            self.treat_prob = sc.smoothinterp(sim.yearvec, self.treat_prob_data.year.values, self.treat_prob_data.treat_prob.values, smoothness=0)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_care_seekers', dtype=int, label="Care seekers", auto_plot=False),
            ss.Result('new_tx0', dtype=int, label="No treatment", auto_plot=False),
            ss.Result('new_tx1', dtype=int, label="1 treatment", auto_plot=False),
            ss.Result('new_tx2', dtype=int, label="2 treatment", auto_plot=False),
            ss.Result('new_tx3', dtype=int, label="3 treatments", auto_plot=False),
        )
        return

    def step(self, uids=None):
        """ Apply syndromic management """
        sim = self.sim
        self.treated_by_uid = None
        ti = self.ti

        # If this intervention has stopped, reset eligibility for all associated treatments
        if sim.now >= self.stop:
            for treatment in self.treatments:
                treatment.eligibility = ss.uids()  # Reset
            return

        if sim.now >= self.start:

            if uids is None:
                uids = self.check_eligibility()
                self.ti_tested[uids] = self.ti

            if len(uids):
                treated_by_uid = np.zeros((len(uids), len(self.treatments)), dtype=bool)

                f_uids = sim.people.female[uids]
                m_uids = sim.people.male[uids]

                # Loop over diseases
                for disease in self.diseases:

                    # Pull out parameters, initialize probabilities
                    sens = self.pars.sens[disease.name]
                    spec = self.pars.spec[disease.name]
                    treatable = disease.treatable[uids]
                    susceptible = disease.susceptible[uids]

                    p_sens = np.full(np.count_nonzero(treatable), np.nan)
                    p_spec = np.full(np.count_nonzero(susceptible), np.nan)

                    inf_f = f_uids[treatable]   # Treatable includes exposed + infected
                    sus_f = f_uids[susceptible]
                    inf_m = m_uids[treatable]
                    sus_m = m_uids[susceptible]

                    p_sens[inf_f] = sens[0]
                    p_spec[sus_f] = spec[0]
                    if disease.name != 'bv':
                        p_sens[inf_m] = sens[1]
                        p_spec[sus_m] = spec[1]

                    # Apply the test
                    self.pars.sens_dist.set(p_sens)
                    self.pars.spec_dist.set(p_spec)
                    true_pos, false_neg = self.pars.sens_dist.split(uids[treatable])
                    true_neg, false_pos = self.pars.spec_dist.split(uids[susceptible])
                    treat_uids = true_pos | false_pos
                    do_treat = np.array([True if u in treat_uids else False for u in uids])

                    # Add to results
                    disease.results['new_true_pos'][ti] += len(true_pos)
                    disease.results['new_false_pos'][ti] += len(false_pos)
                    disease.results['new_true_neg'][ti] += len(true_neg)
                    disease.results['new_false_neg'][ti] += len(false_neg)

                    tx = self.disease_treatment_map[disease.name]
                    if tx is not None:
                        tx_ind = self.treatments.index(tx)
                        treated_by_uid[:, tx_ind] = treated_by_uid[:, tx_ind] | do_treat
                        tx.eligibility = tx.eligibility | treat_uids

                # Update states: time referred to treatment for anyone referred
                referred_uids = uids[treated_by_uid.any(axis=1)]
                dismissed_uids = uids.remove(referred_uids)

                # For females who test negative for NG/CT/TV, refer some to BV treatment
                neg_f = dismissed_uids[self.sim.people.female[dismissed_uids]]
                mtnz = self.pars.p_mtnz.filter(neg_f)
                for neg_tx in self.negative_treatments:
                    neg_tx.eligibility = neg_tx.eligibility | mtnz

                # Update referred and dismissed
                referred_uids = referred_uids | mtnz
                dismissed_uids = uids.remove(referred_uids)
                self.ti_referred[referred_uids] = self.ti
                self.ti_dismissed[dismissed_uids] = self.ti

                # Calculate number of treatments
                n_treatments = treated_by_uid.sum(axis=1)
                mtnz_inds = np.searchsorted(uids, mtnz)
                n_treatments[mtnz_inds] += 1
                self.n_treatments = n_treatments

            return

    def update_results(self):
        super().update_results()
        ti = self.ti
        just_tested = self.ti_tested == ti
        self.results['new_care_seekers'][ti] += count(just_tested)

        # Record the number of people who received 0-3 treatments
        self.results['new_tx0'][ti] += count(self.n_treatments == 0)
        self.results['new_tx1'][ti] += count(self.n_treatments == 1)
        self.results['new_tx2'][ti] += count(self.n_treatments == 2)
        self.results['new_tx3'][ti] += count(self.n_treatments == 3)

        return


class SyndromicMgmt(STITest):
    """
    Syndromic management of vaginal/urethral discharge syndromes.

    When a symptomatic patient presents, they are probabilistically assigned
    one of four treatment outcomes based on whether they have a cervical
    infection (by default NG/CT) and their sex:

        all3  — treated for cervical pathogens + TV/BV
        ngct  — treated for cervical pathogens only
        mtnz  — treated with metronidazole only (TV/BV)
        none  — dismissed without treatment

    Treatment probabilities for each outcome are given separately for
    females with cervical infection (``tx_mix_cerv``), females without
    (``tx_mix_noncerv``), and males (cervical column of ``tx_mix_cerv``
    is reused for UDS).  Each dict maps outcome name → [female_prob, male_prob].

    The ``outcome_treatment_map`` links each outcome to the treatment
    module(s) that should be triggered.  Diagnostic accuracy (TP/FP/FN/TN)
    is recorded on any disease module that exposes those result keys.

    Args:
        treatments (list):           treatment module instances
        diseases (list):             disease modules to track accuracy for
        cervical_diseases (list):    disease modules whose symptomatic flag
                                     indicates cervical infection; defaults to
                                     any disease named 'ng' or 'ct' in ``diseases``
        outcome_treatment_map (dict): maps outcome name → list of treatment modules;
                                     auto-constructed from ``treatments`` if omitted
        treat_prob_data:             unused; reserved for future time-varying tx prob
        years, start, stop, eligibility, name, label: passed to STITest

    Examples::

        syndromic = sti.SyndromicMgmt(
            diseases=[ng, ct, tv, bv],
            treatments=[ng_tx, ct_tx, metronidazole],
            eligibility=seeking_care_vds,
        )
    """

    def __init__(self, pars=None, treatments=None, diseases=None,
                 cervical_diseases=None, outcome_treatment_map=None,
                 treat_prob_data=None, years=None, start=None, stop=None,
                 eligibility=None, name=None, label=None, **kwargs):
        super().__init__(years=years, start=start, stop=stop, eligibility=eligibility,
                         name=name, label=label)
        self.define_pars(
            tx_mix_cerv=dict(
                all3=[0.50, 0.05],
                ngct=[0.20, 0.80],
                mtnz=[0.20, 0.00],
                none=[0.10, 0.15],
            ),
            tx_mix_noncerv=dict(
                all3=[0.40, 0.05],
                ngct=[0.10, 0.80],
                mtnz=[0.20, 0.00],
                none=[0.30, 0.15],
            ),
            tx_cerv_f=ss.choice(a=4),
            tx_cerv_m=ss.choice(a=4),
            tx_noncerv_f=ss.choice(a=4),
            tx_noncerv_m=ss.choice(a=4),
            dt_scale=False,
        )
        self.update_pars(pars, **kwargs)
        self.fvals_cerv    = [v[0] for v in self.pars.tx_mix_cerv.values()]
        self.mvals_cerv    = [v[1] for v in self.pars.tx_mix_cerv.values()]
        self.fvals_noncerv = [v[0] for v in self.pars.tx_mix_noncerv.values()]
        self.mvals_noncerv = [v[1] for v in self.pars.tx_mix_noncerv.values()]

        self.treatments = sc.tolist(treatments)
        self.diseases = sc.tolist(diseases)
        self._cervical_diseases = sc.tolist(cervical_diseases)  # resolved in init_pre

        if outcome_treatment_map is None:
            outcome_treatment_map = dict(
                all3=self.treatments,
                ngct=[self.treatments[0], self.treatments[1]],
                mtnz=[self.treatments[1]],
                none=[],
            )
        self.outcome_treatment_map = outcome_treatment_map

        self.define_states(
            ss.FloatArr('ti_referred'),
            ss.FloatArr('ti_dismissed'),
        )
        self.treat_prob_data = treat_prob_data
        self.treated_by_uid = None
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.pars.tx_cerv_f.set(p=self.fvals_cerv)
        self.pars.tx_cerv_m.set(p=self.mvals_cerv)
        self.pars.tx_noncerv_f.set(p=self.fvals_noncerv)
        self.pars.tx_noncerv_m.set(p=self.mvals_noncerv)
        # Resolve cervical diseases: use explicit list, else find 'ng'/'ct' in diseases
        if not self._cervical_diseases:
            self._cervical_diseases = [d for d in self.diseases if d.name in ('ng', 'ct')]
        return

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        for sk in ['', 'f', 'm']:
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' - {sk.upper()}'
            results += [
                ss.Result('new_care_seekers' + skk, dtype=int, label='Care seekers' + skl),
                ss.Result('new_tx0' + skk, dtype=int, label='No treatment' + skl),
                ss.Result('new_tx1' + skk, dtype=int, label='1 treatment' + skl),
                ss.Result('new_tx2' + skk, dtype=int, label='2 treatments' + skl),
                ss.Result('new_tx3' + skk, dtype=int, label='3 treatments' + skl),
            ]
        self.define_results(*results)
        return

    def step(self, uids=None):
        sim = self.sim
        ppl = sim.people
        self.treated_by_uid = None

        if sim.now >= self.stop:
            for treatment in self.treatments:
                treatment.eligibility = ss.uids()
            return

        if sim.now >= self.start:
            if uids is None:
                uids = self.check_eligibility()
                self.ti_tested[uids] = self.ti

            if len(uids):
                f_uids = uids[ppl.female[uids]]
                m_uids = uids[ppl.male[uids]]

                # Cervical infection: symptomatic in any cervical disease
                if self._cervical_diseases:
                    is_cerv = self._cervical_diseases[0].symptomatic.copy()
                    for d in self._cervical_diseases[1:]:
                        is_cerv = is_cerv | d.symptomatic
                else:
                    is_cerv = ss.BoolArr('is_cerv')
                    is_cerv.initialize(sim.people)

                f_cerv_uids    = f_uids[is_cerv[f_uids]]
                f_noncerv_uids = f_uids[~is_cerv[f_uids]]

                ofc  = self.pars.tx_cerv_f.rvs(f_cerv_uids)
                ofnc = self.pars.tx_noncerv_f.rvs(f_noncerv_uids)
                om   = self.pars.tx_cerv_m.rvs(m_uids)

                outcomes = dict(
                    all3=f_cerv_uids[ofc == 0] | f_noncerv_uids[ofnc == 0] | m_uids[om == 0],
                    ngct=f_cerv_uids[ofc == 1] | f_noncerv_uids[ofnc == 1] | m_uids[om == 1],
                    mtnz=f_cerv_uids[ofc == 2] | f_noncerv_uids[ofnc == 2] | m_uids[om == 2],
                    none=f_cerv_uids[ofc == 3] | f_noncerv_uids[ofnc == 3] | m_uids[om == 3],
                )

                # Diagnostic accuracy: record TP/FP/FN/TN on each disease module
                # all3 and ngct both treat cervical diseases; mtnz treats the rest
                for disease in self.diseases:
                    if not hasattr(disease, 'sex_keys'):
                        continue
                    in_all3_or_ngct = outcomes['all3'] | outcomes['ngct']
                    in_mtnz         = outcomes['mtnz']
                    is_cerv_disease = disease in self._cervical_diseases
                    treat_outcomes  = in_all3_or_ngct if is_cerv_disease else in_mtnz
                    miss_outcomes   = in_mtnz         if is_cerv_disease else in_all3_or_ngct
                    for pkey, pattr in disease.sex_keys.items():
                        skk = '' if pkey == '' else f'_{pkey}'
                        disease.results[f'new_true_pos{skk}'][self.ti]  += len(outcomes['all3'] & disease.treatable & ppl[pattr])
                        disease.results[f'new_false_pos{skk}'][self.ti] += len(outcomes['all3'] & disease.susceptible & ppl[pattr])
                        disease.results[f'new_true_neg{skk}'][self.ti]  += len(outcomes['none'] & disease.susceptible & ppl[pattr])
                        disease.results[f'new_false_neg{skk}'][self.ti] += len(outcomes['none'] & disease.treatable & ppl[pattr])
                        disease.results[f'new_true_pos{skk}'][self.ti]  += len(treat_outcomes & disease.treatable & ppl[pattr])
                        disease.results[f'new_false_pos{skk}'][self.ti] += len(treat_outcomes & disease.susceptible & ppl[pattr])
                        disease.results[f'new_false_neg{skk}'][self.ti] += len(miss_outcomes & disease.treatable & ppl[pattr])
                        disease.results[f'new_true_neg{skk}'][self.ti]  += len(miss_outcomes & disease.susceptible & ppl[pattr])

                # Route to treatments
                for outcome, txs in self.outcome_treatment_map.items():
                    for tx in txs:
                        tx.eligibility = tx.eligibility | outcomes[outcome]

                referred_uids  = outcomes['all3'] | outcomes['ngct'] | outcomes['mtnz']
                dismissed_uids = outcomes['none']
                self.ti_referred[referred_uids]  = self.ti
                self.ti_dismissed[dismissed_uids] = self.ti
                self.treated_by_uid = outcomes

            self.store_results()
        return

    def store_results(self):
        ti  = self.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti
        self.results['new_care_seekers'][ti]   += count(just_tested)
        self.results['new_care_seekers_f'][ti] += count(just_tested & ppl.female)
        self.results['new_care_seekers_m'][ti] += count(just_tested & ppl.male)

        sexdict = {'': 'alive', 'f': 'female', 'm': 'male'}
        if self.treated_by_uid is not None:
            for sk, sl in sexdict.items():
                skk = '' if sk == '' else f'_{sk}'
                self.results['new_tx0' + skk][ti] += count(ppl[sl][self.treated_by_uid['none']])
                self.results['new_tx1' + skk][ti] += count(ppl[sl][self.treated_by_uid['mtnz']])
                self.results['new_tx2' + skk][ti] += count(ppl[sl][self.treated_by_uid['ngct']])
                self.results['new_tx3' + skk][ti] += count(ppl[sl][self.treated_by_uid['all3']])
        return


class ANCTest(ss.Intervention):
    """
    Multi-disease ANC visit testing, scheduled once per pregnancy.

    When a woman becomes pregnant, schedules one ANC visit at a random
    gestational timepoint (months 1–7). At that visit she is tested for all
    supported diseases found in the sim. ``visit_prob`` reflects ANC attendance
    rather than test acceptance — once she attends, all active disease tests are
    applied.

    Auto-detected diseases (registry, Option A):
        hiv      — sets ``hiv.diagnosed`` / ``hiv.ti_diagnosed``; optionally schedules
                   ``InfantHIVTest`` on the unborn via ``maternalnet``
        syphilis — routes positives to ``disease_treatment_map['syphilis']``; optionally
                   schedules ``NewbornSyphTest`` on the unborn via ``maternalnet``
    Additional diseases in ``disease_names`` are tested and routed via
    ``disease_treatment_map`` with configurable sensitivity.

    Newborn/infant tests must be added to ``sim.interventions`` by the caller and
    passed via ``newborn_tests`` so they go through normal Starsim initialisation.
    ``ANCTest`` only handles scheduling them (sets their ``ti_scheduled`` state).

    Args:
        visit_prob (float):           probability of attending ANC (default 0.85)
        disease_names (list):         diseases to test for; default: auto-detect
                                      'hiv' and 'syphilis' from sim
        test_sensitivity (dict):      per-disease sensitivity, e.g. ``{'hiv': 0.99}``;
                                      default 1.0 for all
        disease_treatment_map (dict): maps disease name → treatment intervention;
                                      required for syphilis and any non-HIV disease
        newborn_tests (dict):         maps disease name → newborn/infant test
                                      intervention already in ``sim.interventions``,
                                      e.g. ``{'hiv': infant_hiv, 'syphilis': newborn_syph}``
        start, stop, name, label:     standard

    Examples::

        # Minimal: HIV + syphilis auto-detected, no newborn tests
        anc = sti.ANCTest(
            visit_prob=0.85,
            disease_treatment_map={'syphilis': syph_tx},
        )

        # With newborn/infant tests
        infant_hiv  = sti.InfantHIVTest(name='infant_hiv')
        newborn_syph = sti.NewbornSyphTest(name='newborn_syph', ...)
        anc = sti.ANCTest(
            visit_prob=0.85,
            disease_treatment_map={'syphilis': syph_tx},
            newborn_tests={'hiv': infant_hiv, 'syphilis': newborn_syph},
        )
        sim = ss.Sim(..., interventions=[anc, art, syph_tx, infant_hiv, newborn_syph])
    """

    # Registry: disease names → (needs_diagnosed_flag, has_newborn_test)
    # HIV is handled via diagnosed state; syphilis via treatment map.
    # Lazy imports for the actual classes live in _schedule_newborn().
    _registry = {
        'hiv':      dict(use_diagnosed=True),
        'syphilis': dict(use_diagnosed=False),
    }

    def __init__(self, visit_prob=0.85, disease_names=None, test_sensitivity=None,
                 disease_treatment_map=None, newborn_tests=None,
                 start=None, stop=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label, **kwargs)
        self._visit_prob_val  = visit_prob
        self._disease_names   = disease_names          # None → auto-detect from registry
        self._test_sensitivity = test_sensitivity or {}
        self._disease_tx_map  = disease_treatment_map or {}
        self._newborn_tests   = newborn_tests or {}
        self._active_diseases = []   # resolved in init_pre
        self._visit_prob_dist = ss.bernoulli(p=visit_prob)
        self._sens_dist       = ss.bernoulli(p=0.5)
        self._visit_timing    = ss.randint(low=1, high=9)   # months into pregnancy
        self.start = start
        self.stop  = stop
        self.define_states(ss.FloatArr('ti_visit', label='Scheduled ANC visit'))
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start is None:
            self.start = sim.t.yearvec[0]
        if self.stop is None:
            self.stop = sim.t.yearvec[-1]

        if self._disease_names is not None:
            self._active_diseases = list(self._disease_names)
        else:
            # Auto-detect: include any registered disease present in the sim
            self._active_diseases = [d for d in self._registry if d in sim.diseases]
        return

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        results += ss.Result('n_attended', dtype=int, label='ANC attendees')
        for dname in self._active_diseases:
            results += ss.Result(f'n_{dname}_positive', dtype=int,
                                 label=f'{dname.upper()} positive at ANC')
        self.define_results(*results)
        return

    def step(self):
        sim = self.sim
        ti  = self.ti

        if sim.now < self.start or sim.now >= self.stop:
            return

        # Schedule a visit for each newly pregnant woman
        if hasattr(sim.demographics, 'pregnancy'):
            newly_preg = (sim.demographics.pregnancy.ti_pregnant == ti).uids
            attendees  = self._visit_prob_dist.filter(newly_preg)
            months_dt  = ss.dur(1, 'month') / sim.t.dt   # timesteps per month
            offsets    = np.round(self._visit_timing.rvs(attendees) * months_dt).astype(int)
            self.ti_visit[attendees] = ti + offsets

        visit_uids = (self.ti_visit == ti).uids
        if not len(visit_uids):
            return

        self.results['n_attended'][ti] = len(visit_uids)

        for dname in self._active_diseases:
            if dname not in sim.diseases:
                continue
            disease  = sim.diseases[dname]
            sens     = self._test_sensitivity.get(dname, 1.0)

            # Apply sensitivity to infected attendees
            infected = visit_uids[disease.infected[visit_uids]]
            if len(infected) and sens < 1.0:
                self._sens_dist.set(p=sens)
                pos_uids = infected[self._sens_dist.rvs(infected)]
            else:
                pos_uids = infected

            if not len(pos_uids):
                continue

            self.results[f'n_{dname}_positive'][ti] += len(pos_uids)
            cfg = self._registry.get(dname, {})

            if cfg.get('use_diagnosed'):
                # HIV-style: set diagnosed flag so ART picks up automatically
                disease.diagnosed[pos_uids]    = True
                disease.ti_diagnosed[pos_uids] = ti
            elif dname in self._disease_tx_map:
                self._disease_tx_map[dname].eligibility = (
                    self._disease_tx_map[dname].eligibility | pos_uids
                )

            # Schedule newborn/infant test if provided and maternalnet exists
            if dname in self._newborn_tests:
                self._schedule_newborn(sim, pos_uids, self._newborn_tests[dname])
        return

    @staticmethod
    def _schedule_newborn(sim, pos_mother_uids, newborn_test):
        """Schedule a newborn/infant test for unborn children of positive mothers."""
        if not (hasattr(sim.networks, 'maternalnet') and
                hasattr(sim.demographics, 'pregnancy')):
            return
        mn = sim.networks.maternalnet
        pos_mask = np.isin(mn.p1, pos_mother_uids)
        unborn   = mn.p2[pos_mask]
        ti_births = sim.demographics.pregnancy.ti_delivery[
            mn.p1[pos_mask]
        ]
        valid = ~np.isnan(ti_births)
        if valid.any():
            newborn_test.schedule(unborn[valid], ti_births[valid].astype(int))
        return


class STITreatment(ss.Intervention):
    """
    Base class for treatment of STI infection.
    The majority of STI treatments will clear infection.

    Args:
        pars:
        disease (str): should match the name of one of the diseases in the simulation
    """
    def __init__(self, name=None, pars=None, diseases=None, eligibility=None, max_capacity=None, years=None, *args, **kwargs):
        super().__init__(*args, name=name)
        self.requires = diseases
        self.define_pars(
            treat_prob=ss.bernoulli(p=1.),
            treat_eff=ss.bernoulli(p=0.9),
            by_sex=True,  # Whether or not to store outcomes by sex
        )
        self.diseases = sc.promotetolist(diseases)
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        if self.eligibility is None:
            self.eligibility = ss.uids()
        self.queue = []
        self.max_capacity = max_capacity
        self.years = years

        # Results by sex - set during init_results
        self.sexkeys = None
        self.sexdict = None

        # States
        self.define_states(
            ss.BoolState('treated'),
            ss.FloatArr('ti_treated'),
        )

        return

    def init_results(self):
        super().init_results()
        if self.pars.by_sex:
            self.sexkeys = ['', 'f', 'm']
            self.sexdict = {'': 'alive', 'f': 'female', 'm': 'male'}
        else:
            self.sexkeys = ['']
            self.sexdict = {'': 'alive'}
        results = sc.autolist()
        for sk in self.sexkeys:
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' - {sk.upper()}'
            results += [
                ss.Result('new_treated'+skk, dtype=int, label="Number treated"+skl, auto_plot=False),
                ss.Result('new_treated_success'+skk, dtype=int, label="Successfully treated"+skl, auto_plot=False),
                ss.Result('new_treated_failure'+skk, dtype=int, label="Treatment failure"+skl, auto_plot=False),
                ss.Result('new_treated_unnecessary'+skk, dtype=int, label="Overtreatment"+skl, auto_plot=False),
            ]
        self.define_results(*results)
        return

    def get_candidates(self, sim):
        """
        Get people who are willing to accept treatment
        """
        accept_uids = ss.uids()
        if callable(self.eligibility):
            eligible_uids = self.check_eligibility()  # Apply eligiblity - uses base class from ss.Intervention
        else:
            eligible_uids = self.eligibility
        if len(eligible_uids):
            accept_uids = self.pars.treat_prob.filter(eligible_uids)
        return accept_uids

    def set_treat_eff(self, uids):
        """ Can be changed by derived classes """
        return

    def administer(self, sim, uids, disease, return_format='dict'):
        """ Administer treatment, keeping track of unnecessarily treated individuals """

        inf = uids & sim.diseases[disease].treatable
        sus = uids & sim.diseases[disease].susceptible

        self.set_treat_eff(inf)
        successful = self.pars.treat_eff.filter(inf)
        unsuccessful = np.setdiff1d(inf, successful)
        unnecessary = sus

        if self.pars.by_sex:
            successful_f = successful[sim.people.female[successful]]
            unsuccessful_f = unsuccessful[sim.people.female[unsuccessful]]
            unnecessary_f = unnecessary[sim.people.female[unnecessary]]
            successful_m = successful[sim.people.male[successful]]
            unsuccessful_m = unsuccessful[sim.people.male[unsuccessful]]
            unnecessary_m = unnecessary[sim.people.male[unnecessary]]

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
            if self.pars.by_sex:
                output = sc.mergedicts(output,
                    {'successful_f': successful_f, 'unsuccessful_f': unsuccessful_f, 'unnecessary_f': unnecessary_f,
                     'successful_m': successful_m, 'unsuccessful_m': unsuccessful_m, 'unnecessary_m': unnecessary_m})
        elif return_format == 'array':
            output = successful

        return output

    def change_states(self, disease, treat_succ):
        """ Change the states of people who are treated """
        self.sim.diseases[disease].clear_infection(treat_succ)
        self.sim.diseases[disease].wipe_dates(treat_succ)
        return

    def step(self):
        """
        Apply treatment. On each timestep, this method will add eligible people who are willing to accept treatment to a
        queue, and then will treat as many people in the queue as there is capacity for.
        """
        sim = self.sim
        self.outcomes = sc.objdict()
        for disease in self.diseases:
            self.outcomes[disease] = sc.objdict(successful=ss.uids(), unsuccessful=ss.uids(), unnecessary=ss.uids())
            if self.pars.by_sex:
                self.outcomes[disease].update(successful_f=ss.uids(), unsuccessful_f=ss.uids(), unnecessary_f=ss.uids(),
                                              successful_m=ss.uids(), unsuccessful_m=ss.uids(), unnecessary_m=ss.uids())

        # Figure out who to treat
        treat_uids = self.get_candidates(sim)
        self.treated[treat_uids] = True
        self.ti_treated[treat_uids] = self.ti

        # Treat people
        if len(treat_uids):
            for disease in self.diseases:
                self.outcomes[disease] = self.administer(sim, treat_uids, disease)

                # Change states
                treat_succ = self.outcomes[disease]['successful']
                if len(treat_succ):
                    self.change_states(disease, treat_succ)

            # If eligibility is a list of UIDS, remove the treated ones. WARNING, FRAGILE!!!
            if isinstance(self.eligibility, ss.uids):
                self.eligibility = self.eligibility.remove(treat_uids)

        successful = ss.uids()
        unsuccessful = sc.dcp(treat_uids)
        unneeded = sc.dcp(treat_uids)

        for disease in self.diseases:
            # Unneeded if it doesn't do anything for any disease
            unneeded = unneeded.remove(self.outcomes[disease]['successful'] | self.outcomes[disease]['unsuccessful'])
            # Successful if it's successful for at least one disease
            successful = successful | self.outcomes[disease]['successful']
            # Unsuccessful if it's unsuccessful for all diseases
            unsuccessful = unsuccessful & self.outcomes[disease]['unsuccessful']

        for sk, sl in self.sexdict.items():
            skk = '' if sk == '' else f'_{sk}'
            self.outcomes['unnecessary'+skk] = unneeded[sim.people[sl][unneeded]]
            self.outcomes['successful'+skk] = successful[sim.people[sl][successful]]
            self.outcomes['unsuccessful'+skk] = unsuccessful[sim.people[sl][unsuccessful]]

        return treat_uids

    def update_results(self):
        ti = self.ti
        just_treated = self.ti_treated == ti

        # Store new treatment results in the disease module results
        for sk, sl in self.sexdict.items():
            skk = '' if sk == '' else f'_{sk}'
            self.results['new_treated_success'+skk][ti] = len(self.outcomes['successful'+skk])
            self.results['new_treated_failure'+skk][ti] = len(self.outcomes['unsuccessful'+skk])
            self.results['new_treated_unnecessary'+skk][ti] = len(self.outcomes['unnecessary'+skk])
            self.results['new_treated'+skk][ti] = count(just_treated & self.sim.people[sl])
            for disease in self.diseases:
                try:
                    self.sim.diseases[disease].results['new_treated_success'+skk][ti] += len(self.outcomes[disease]['successful'+skk])
                    self.sim.diseases[disease].results['new_treated_failure'+skk][ti] += len(self.outcomes[disease]['unsuccessful'+skk])
                    self.sim.diseases[disease].results['new_treated_unnecessary'+skk][ti] += len(self.outcomes[disease]['unnecessary'+skk])
                    self.sim.diseases[disease].results['new_treated'+skk][ti] += count(just_treated & self.sim.people[sl])
                except:
                    print('error')

        # Debugging
        if self.name in ['ng_tx', 'ct_tx']:
            for oc in ['new_treated_success', 'new_treated_failure', 'new_treated_unnecessary']:
                if self.sim.diseases[disease].results[oc][ti] != self.results[oc][ti]:
                    errormsg = 'Inconsistent treatment outcomes'
                    raise ValueError(errormsg)

        return


class PartnerNotification(ss.Intervention):

    def __init__(self, eligibility, test, test_prob=0.5, **kwargs):
        """
        :param disease: The disease module from which to draw the transmission tree used to find contacts
        :param eligible: A function `f(sim)` that returns the UIDs/BoolArr of people to trace (typically people who just tested positive)
        :param test: The testing intervention to use when testing identified contacts
        :param test_prob: The probability of a contact being identified and accepting a test
        :param kwargs: Other arguments passed to ``ss.Intervention``
        """
        super().__init__(**kwargs)
        self.eligibility = eligibility
        self.test = test
        self.test_prob = ss.bernoulli(test_prob)

    def identify_contacts(self, uids):
        # Return UIDs of people that have been identified as contacts and should be notified

        if len(uids) == 0:
            return ss.uids()

        # Find current contacts
        nw = self.sim.networks.structuredsexual
        contacts = nw.find_contacts(uids, as_array=False)  # Current contacts

        # Find historical contacts
        log = self.sim.diseases[self.disease].log
        for source, _, network in log.in_edges(uids, data="network"):
            if network == 'structuredsexual':
                contacts.add(source)  # Add the infecting agents

        for _, target, network in log.out_edges(uids, data="network"):
            if network == 'structuredsexual':
                contacts.add(target)  # Add infected agents

        # Filter by test_prob and return UIDs
        return self.test_prob.filter(ss.uids(contacts))

    def notify(self, uids):
        # Schedule a test for identified contacts at the next timestep (this also ensures that contacts tracing will take place for partners that test positive)
        # Could include a parameter here for acceptance of testing (if separating out probabilities of notification and testing)
        return self.test.schedule(uids, self.ti+1)

    def step(self):
        uids = self.eligibility(self.sim)
        uids = self.identify_contacts(uids)
        return self.notify(uids)


class ProductMix(ss.Product):
    """
    Generic class for algorithms that determine which product a person will receive
    Uses ss.choice() sampling, which is slower than bernoulli, when there are more than two options
    The test that agents are given does NOT depend on their underlying health state.
    """

    def __init__(self, df, excl_cols=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.df = df
        self.result_list = df.products.unique()
        self.product_dist = ss.choice(a=self.result_list)
        self.product_mix = None  # Set during initialization

        # Pull out data from dataframe
        if excl_cols is None:
            excl_cols = 'products'

        self.data_years = np.array([float(c) for c in df.columns if c not in excl_cols])
        data_years_str = [c for c in df.columns if c not in excl_cols]
        y = df[data_years_str].values
        self.f_out = interp1d(self.data_years, y, axis=1, fill_value=(y[:, 0], y[:, -1]), bounds_error=False, kind='zero')
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.product_mix = self.f_out(self.t.yearvec)

    def administer(self, sim, uids):
        """
        Apply a testing algorithm
        """
        outcomes = {r: ss.uids() for r in self.result_list}
        self.product_dist.set(p=self.product_mix[:, self.ti])
        this_result = self.product_dist.rvs(uids)
        if len(this_result):
            for res in self.result_list:
                outcomes[res] = uids[this_result == res]

        return outcomes


