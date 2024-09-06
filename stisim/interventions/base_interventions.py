"""
Define interventions for STIsim
"""

import starsim as ss
import numpy as np
import sciris as sc
from scipy.interpolate import interp1d


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


# %% Base classes
__all__ = ["STIDx", "STITest", "SyndromicMgmt", "STITreatment", "PartnerNotification", "ProductMix"]


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

        # Store the probability of returning each possible result
        self.result_dist = dict()
        for state in self.health_states:
            thisdf = self.df.loc[self.df.state == state]
            self.result_dist[state] = ss.bernoulli(thisdf.p_positive.values)
        return

    def administer(self, sim, uids):
        """
        Administer a testing product.
        """
        outcomes = {r: ss.uids() for r in self.result_list}
        for state in self.health_states:
            this_state = getattr(sim.diseases[self.disease], state)
            true_uids = this_state.uids  # Find people for which this state is true
            if len(true_uids):
                these_uids = true_uids.intersect(uids)  # Find intersection of people in this state and supplied UIDs
                pos, neg = self.result_dist[state].split(these_uids)
                if len(pos): outcomes['positive'] = outcomes['positive'] | pos
                if len(neg): outcomes['negative'] = outcomes['negative'] | neg
        return outcomes


class STITest(ss.Intervention):
    """
    Base class for STI testing
    """

    def __init__(self, pars=None, test_prob_data=None, years=None, start=None, end=None, eligibility=None, product=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.default_pars(
            rel_test=1,
        )
        self.update_pars(pars, **kwargs)

        # Years
        if years is not None and start is not None:
            errormsg = 'Provide either years or start_year, not both.'
            raise ValueError(errormsg)
        self.years = years
        self.start = start
        self.end = end
        if self.start is None:
            if self.years is not None:
                self.start = self.years[0]
        if self.end is None:
            if self.years is not None:
                self.end = self.years[-1]

        # Testing eligibility and uptake
        self.eligibility = eligibility
        self.test_prob_data = test_prob_data
        self.test_prob = ss.bernoulli(self.make_test_prob_fn)

        # States
        self.add_states(
            ss.BoolArr('tested'),
            ss.BoolArr('diagnosed'),
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
            self.outcomes = {outcome: ss.FloatArr(f'ti_{outcome}') for outcome in self.product.result_list}
            self.last_outcomes = {outcome: ss.uids() for outcome in self.product.result_list}

        return

    @property
    def states(self):
        return super().states + list(self.outcomes.values())

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start is None: self.start = sim.pars.start
        if self.end is None: self.end = sim.pars.end
        self.init_results()
        return

    def init_results(self):
        npts = self.sim.npts
        results = [
            ss.Result(self.name, 'new_diagnoses', npts, dtype=float, scale=True, label="New diagnoses"),
            ss.Result(self.name, 'new_tests', npts, dtype=int, scale=True, label="New tests"),
        ]
        self.results += results

        return

    @staticmethod
    def make_test_prob_fn(self, sim, uids):
        """ Testing probabilites over time """

        if sc.isnumber(self.test_prob_data):
            test_prob = self.test_prob_data

        elif sc.checktype(self.test_prob_data, 'arraylike'):
            year_ind = sc.findnearest(self.years, sim.year)
            test_prob = self.test_prob_data[year_ind]
        else:
            errormsg = 'Format of test_prob_data must be float or array.'
            raise ValueError(errormsg)

        # Scale and validate
        test_prob = test_prob * self.pars.rel_test * sim.dt
        test_prob = np.clip(test_prob, a_min=0, a_max=1)

        return test_prob

    def get_testers(self, sim):
        """
        Find who tests by applying eligibility and coverage/uptake
        """
        # Select UIDs for testing based on eligibility and test_prob
        accept_uids = ss.uids()
        eligible_uids = self.check_eligibility(sim)  # Apply eligiblity
        if len(eligible_uids):
            accept_uids = self.test_prob.filter(eligible_uids)
        scheduled_uids = (self.ti_scheduled == sim.ti).uids  # Add on scheduled tests
        return accept_uids | scheduled_uids

    def schedule(self, uids, ti):
        # Schedule a future test for specified UIDs at specified time indices
        self.ti_scheduled[uids] = ti

    def apply(self, sim, uids=None):
        """ Apply the testing intervention """
        outcomes = {outcome: ss.uids() for outcome in self.product.result_list}
        self.last_outcomes = outcomes

        # Apply if within the start years
        if (sim.year >= self.start) & (sim.year < self.end):

            if uids is None:
                uids = self.get_testers(sim)
                self.ti_tested[uids] = sim.ti

            if len(uids):
                outcomes = self.product.administer(sim, uids)
                self.last_outcomes = outcomes
                for k in self.product.result_list:
                    self.outcomes[k][self.last_outcomes[k]] = sim.ti
            else:
                return outcomes

            # Set time of diagnosis or dismissal
            self.ti_negative[outcomes['negative']] = sim.ti
            self.ti_positive[outcomes['positive']] = sim.ti

            # Update results
            self.update_results()

        return outcomes

    def update_results(self):
        # Store results
        ti = self.sim.ti
        self.results['new_diagnoses'][ti] += count(self.ti_positive == ti)
        self.results['new_tests'][ti] += count(self.ti_tested == ti)


class SyndromicMgmt(STITest):
    """
    Base class for syndromic management.
    Unlike other test classes, this doesn't return positive/negative outcomes, since
    syndromic management doesn't involve reaching a diagnosis.
    Rather, the testing intervention itself contains a linked treatment intervention.
    """

    def __init__(self, pars=None, treatments=None, diseases=None, treat_prob_data=None, years=None, start=None, end=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(years=years, start=start, end=end, eligibility=eligibility, name=name, label=label)
        self.default_pars(
            p_treat=ss.bernoulli(p=0.9),  # The probability that a person with discharge will be offered treatment
        )
        self.update_pars(pars, **kwargs)
        self.treatments = sc.tolist(treatments)
        self.diseases = diseases
        self.add_states(
            ss.FloatArr('ti_referred'),
            ss.FloatArr('ti_dismissed'),
        )
        self.treat_prob_data = treat_prob_data
        self.treat_prob = None
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.treat_prob_data is not None:
            self.treat_prob = np.interp(sim.yearvec, self.treat_prob_data.year.values, self.treat_prob_data.treat_prob.values)
        return

    def init_results(self):
        super().init_results()
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'care_seekers', npts, dtype=int, scale=True, label="Care seekers"),
        ]
        for disease in self.diseases:
            results = [
                ss.Result(disease.name, 'new_false_neg', npts, dtype=float, scale=True, label="False negatives"),
                ss.Result(disease.name, 'new_true_neg', npts, dtype=float, scale=True, label="True negatives"),
                ss.Result(disease.name, 'new_false_pos', npts, dtype=float, scale=True, label="False positive"),
                ss.Result(disease.name, 'new_true_pos', npts, dtype=float, scale=True, label="True positives"),
            ]
            disease.results += results
        return

    def apply(self, sim, uids=None):
        """ Apply syndromic management """
        if (sim.year >= self.start) & (sim.year < self.end):

            if uids is None:
                uids = self.check_eligibility(sim)
                self.ti_tested[uids] = sim.ti

            if len(uids):
                if self.treat_prob is not None:
                    self.pars.p_treat.set(self.treat_prob[sim.ti])
                treat_uids, dismiss_uids = self.pars.p_treat.split(uids)
                for treatment in self.treatments:
                    treatment.eligibility = treat_uids
                self.ti_referred[treat_uids] = sim.ti
                self.ti_dismissed[dismiss_uids] = sim.ti

            # Update results
            self.update_results()

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        treat_uids = self.ti_referred == ti
        dismiss_uids = self.ti_dismissed == ti
        just_tested = self.ti_tested == ti
        self.results['care_seekers'][ti] += count(just_tested)
        for disease in self.diseases:
            disease.results['new_true_pos'][ti] += count(disease.treatable[treat_uids])
            disease.results['new_false_pos'][ti] += count(disease.susceptible[treat_uids])
            disease.results['new_true_neg'][ti] += count(disease.susceptible[dismiss_uids])
            disease.results['new_false_pos'][ti] += count(disease.treatable[dismiss_uids])

        return


class STITreatment(ss.Intervention):
    """
    Base class for treatment of STI infection.
    The majority of STI treatments will clear infection.
    Args:
        pars:
        disease (str): should match the name of one of the diseases in the simulation
    """
    def __init__(self, pars=None, disease=None, eligibility=None, max_capacity=None, years=None, *args, **kwargs):
        super().__init__(*args)
        self.requires = disease
        self.default_pars(
            treat_prob=ss.bernoulli(p=0.9),
            treat_eff=ss.bernoulli(p=0.95),
        )
        self.disease = disease
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        if self.eligibility is None:
            self.eligibility = ss.uids()
        self.queue = []
        self.max_capacity = max_capacity
        self.years = years

        # States
        self.add_states(
            ss.BoolArr('treated'),
            ss.FloatArr('ti_treated'),
        )

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.init_results()
        return

    def init_results(self):
        results = [
            ss.Result(self.disease, 'new_treated', self.sim.npts, dtype=int, scale=True, label="Number treated"),
            ss.Result(self.disease, 'new_treated_success', self.sim.npts, dtype=int, scale=True, label="Successfully treated"),
            ss.Result(self.disease, 'new_treated_failure', self.sim.npts, dtype=int, scale=True, label="Treatment failure"),
            ss.Result(self.disease, 'new_treated_unnecessary', self.sim.npts, dtype=int, scale=True, label="Overtreatment"),
        ]
        # self.sim.diseases[self.disease].results += results
        self.results += results
        return

    def add_to_queue(self, sim):
        """
        Add people who are willing to accept treatment to the queue
        """
        accept_uids = ss.uids()
        if callable(self.eligibility):
            eligible_uids = self.check_eligibility(sim)  # Apply eligiblity - uses base class from ss.Intervention
        else:
            eligible_uids = self.eligibility
        if len(eligible_uids):
            accept_uids = self.pars.treat_prob.filter(eligible_uids)
        if len(accept_uids): self.queue += accept_uids.tolist()
        return

    def get_candidates(self, sim):
        """
        Get the indices of people who are candidates for treatment
        """
        treat_candidates = np.array([], dtype=int)

        if len(self.queue):

            if self.max_capacity is None:
                treat_candidates = self.queue[:]

            else:
                if sc.isnumber(self.max_capacity):
                    max_capacity = self.max_capacity
                elif sc.checktype(self.max_capacity, 'arraylike'):
                    year_ind = sc.findnearest(self.years, sim.year)
                    max_capacity = self.max_capacity[year_ind]

                if max_capacity > len(self.queue):
                    treat_candidates = self.queue[:]
                else:
                    treat_candidates = self.queue[:self.max_capacity]

        return ss.uids(treat_candidates)

    def set_treat_eff(self, uids):
        """ Can be changed by derived classes """
        return

    def administer(self, sim, uids, return_format='dict'):
        """ Administer treatment, keeping track of unnecessarily treated individuals """

        inf = sim.diseases[self.disease].treatable
        sus = sim.diseases[self.disease].susceptible
        inf_uids = uids[inf[uids]]
        sus_uids = uids[sus[uids]]

        self.set_treat_eff(inf_uids)
        successful = self.pars.treat_eff.filter(inf_uids)
        unsuccessful = np.setdiff1d(inf_uids, successful)
        unnecessary = sus_uids

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
        elif return_format == 'array':
            output = successful

        return output

    def change_states(self, sim, treat_succ):
        """ Change the states of people who are treated """
        self.sim.diseases[self.disease].clear_infection(treat_succ)
        self.sim.diseases[self.disease].wipe_dates(treat_succ)
        return

    def apply(self, sim):
        """
        Apply treatment. On each timestep, this method will add eligible people who are willing to accept treatment to a
        queue, and then will treat as many people in the queue as there is capacity for.
        """
        self.outcomes = sc.objdict(successful=[], unsuccessful=[], unnecessary=[])

        # Figure out who to treat
        self.add_to_queue(sim)
        treat_uids = self.get_candidates(sim)
        self.treated[treat_uids] = True
        self.ti_treated[treat_uids] = sim.ti

        # Treat people
        if len(treat_uids):
            self.outcomes = self.administer(sim, uids=treat_uids)
        self.queue = [e for e in self.queue if e not in treat_uids]  # Recreate the queue, removing treated people

        # Change states
        treat_succ = self.outcomes['successful']
        if len(treat_succ):
            self.change_states(sim, treat_succ)

        # Update results
        self.update_results()

        return treat_uids

    def update_results(self):
        ti = self.sim.ti
        treat_uids = (self.ti_treated == ti).uids

        # Store new treatment results both in this intervention and in the gonorrhea disease module results
        self.results['new_treated_success'][ti] = len(self.outcomes['successful'])
        self.results['new_treated_failure'][ti] = len(self.outcomes['unsuccessful'])
        self.results['new_treated_unnecessary'][ti] = len(self.outcomes['unnecessary'])
        self.results['new_treated'][ti] = len(treat_uids)
        self.sim.diseases[self.disease].results['new_treated_success'][ti] += len(self.outcomes['successful'])
        self.sim.diseases[self.disease].results['new_treated_failure'][ti] += len(self.outcomes['unsuccessful'])
        self.sim.diseases[self.disease].results['new_treated_unnecessary'][ti] += len(self.outcomes['unnecessary'])
        self.sim.diseases[self.disease].results['new_treated'][ti] += len(treat_uids)

        return


class PartnerNotification(ss.Intervention):

    def __init__(self, disease, eligible, test, test_prob=0.5, **kwargs):
        """

        :param disease: The disease module from which to draw the transmission tree used to find contacts
        :param eligible: A function `f(sim)` that returns the UIDs/BoolArr of people to trace (typically people who just tested positive)
        :param test: The testing intervention to use when testing identified contacts
        :param test_prob: The probability of a contact being identified and accepting a test
        :param kwargs: Other arguments passed to ``ss.Intervention``
        """
        super().__init__(**kwargs)
        self.disease = disease
        self.eligible = eligible
        self.test = test
        self.test_prob = ss.bernoulli(test_prob)

    def identify_contacts(self, sim, uids):
        # Return UIDs of people that have been identified as contacts and should be notified

        if len(uids) == 0:
            return ss.uids()

        # Find current contacts
        contacts = sim.networks['structuredsexual'].find_contacts(uids, as_array=False) # Current contacts

        # Find historical contacts
        log = sim.diseases[self.disease].log
        for source, _, network in log.in_edges(uids, data="network"):
            if network == 'structuredsexual':
                contacts.add(source) # Add the infecting agents

        for _, target, network in log.out_edges(uids, data="network"):
            if network == 'structuredsexual':
                contacts.add(target)  # Add infected agents

        # Filter by test_prob and return UIDs
        return self.test_prob.filter(ss.uids(contacts))


    def notify(self, sim, uids):
        # Schedule a test for identified contacts at the next timestep (this also ensures that contacts tracing will take place for partners that test positive)
        # Could include a parameter here for acceptance of testing (if separating out probabilities of notification and testing)
        # print(f'Scheduling {len(uids)} tests (of whom {(~sim.diseases.syphilis.susceptible[uids]).sum()} will be positive)')
        return self.test.schedule(uids, sim.ti+1)

    def apply(self, sim):
        uids = self.eligible(sim)
        uids = self.identify_contacts(sim, uids)
        return self.notify(sim, uids)


class ProductMix(ss.Product):
    """
    Generic class for algorithms that determine which product a person will receive
    Uses ss.choice() sampling, which is slower than bernoulli, when there are more than two options
    The test that agents are given does NOT depend on their underlying health state.
    """

    def __init__(self, df, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.df = df
        self.result_list = df.products.unique()
        self.product_dist = ss.choice(a=self.result_list)
        self.product_mix = None  # Set during initialization

        # Pull out data from dataframe
        self.data_years = np.array([float(c) for c in df.columns if c != 'products'])
        data_years_str = [c for c in df.columns if c != 'products']
        y = df[data_years_str].values
        self.f_out = interp1d(self.data_years, y, axis=1, fill_value=(y[:, 0], y[:, -1]), bounds_error=False)
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.product_mix = self.f_out(sim.yearvec)

    def administer(self, sim, uids):
        """
        Apply a testing algorithm
        """
        outcomes = {r: ss.uids() for r in self.result_list}
        self.product_dist.set(p=self.product_mix[:, sim.ti])
        this_result = self.product_dist.rvs(uids)
        if len(this_result):
            for res in self.result_list:
                outcomes[res] = uids[this_result == res]

        return outcomes

