"""
Define interventions for STIsim
"""

import starsim as ss
import numpy as np
import sciris as sc
import pandas as pd
from scipy.interpolate import interp1d


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


# %% Base classes
__all__ = ["STIDx", "STITest", "SymptomaticTesting", "STITreatment", "PartnerNotification", "ProductMix"]


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
    Base class for STI testing
    """

    def __init__(self, pars=None, test_prob_data=None, years=None, start=None, stop=None, eligibility=None, product=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            dt='month',
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
            ss.Result('new_diagnoses', dtype=int, label="New diagnoses"),
            ss.Result('new_tests', dtype=int, label="New tests"),
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
        if self.pars.dt_scale: test_prob *= self.t.dt_year
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
            self.treat_prob = np.interp(sim.yearvec, self.treat_prob_data.year.values, self.treat_prob_data.treat_prob.values)
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_care_seekers', dtype=int, label="Care seekers"),
            ss.Result('new_tx0', dtype=int, label="No treatment"),
            ss.Result('new_tx1', dtype=int, label="1 treatment"),
            ss.Result('new_tx2', dtype=int, label="2 treatment"),
            ss.Result('new_tx3', dtype=int, label="3 treatments"),
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
            dt='month',
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
                ss.Result('new_treated'+skk, dtype=int, label="Number treated"+skl),
                ss.Result('new_treated_success'+skk, dtype=int, label="Successfully treated"+skl),
                ss.Result('new_treated_failure'+skk, dtype=int, label="Treatment failure"+skl),
                ss.Result('new_treated_unnecessary'+skk, dtype=int, label="Overtreatment"+skl),
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


