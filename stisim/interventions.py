"""
Define interventions for STIsim
"""

import starsim as ss
import numpy as np
import sciris as sc


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


# %% Base classes
__all__ = ["STIDx", "STITest", "SyndromicMgmt", "STITreatment", "PartnerNotification"]


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
            ss.FloatArr('ti_diagnosed'),
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

            # Set time of diagnosis
            self.ti_diagnosed[outcomes['positive']] = sim.ti

            # Update results
            self.update_results()

        return outcomes

    def update_results(self):
        # Store results
        ti = self.sim.ti
        self.results['new_diagnoses'][ti] += count(self.ti_diagnosed == ti)
        self.results['new_tests'][ti] += count(self.ti_tested == ti)


class SyndromicMgmt(STITest):
    """
    Base class for syndromic management.
    Unlike other test classes, this doesn't return positive/negative outcomes, since
    syndromic management doesn't involve reaching a diagnosis.
    Rather, the testing intervention itself contains a linked treatment intervention.
    """

    def __init__(self, pars=None, treatments=None, diseases=None, test_prob_data=None, years=None, start=None, end=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, end=end, eligibility=eligibility, name=name, label=label)
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
        self.sim.diseases[self.disease].results += results
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

    @staticmethod
    def change_states(sim, treat_succ):
        """ Change the states of people who are treated. Must be implemented by derived classes """
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


# %% HIV classes
__all__ += ["HIVDx", "HIVTest", "ART", "VMMC"]


class HIVDx(ss.Product):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.result_list = ['positive', 'negative']

    def administer(self, sim, uids):
        outcomes = {r: ss.uids() for r in self.result_list}
        outcomes['positive'] = sim.diseases.hiv.infected.uids.intersect(uids)
        outcomes['negative'] = sim.diseases.hiv.susceptible.uids.intersect(uids)
        return outcomes


class HIVTest(STITest):
    """
    Base class for HIV testing
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None, eligibility=None, name=None, label=None, **kwargs):
        if product is None: product = HIVDx()
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years, start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed

    def apply(self, sim, uids=None):
        outcomes = super().apply(sim, uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = sim.ti
        return outcomes

class ART(ss.Intervention):
    """
    ART-treatment intervention by Robyn Stuart, Daniel Klein and Cliff Kerr, edited by Alina Muellenmeister
    """

    def __init__(self, pars=None, coverage_data=None, start_year=None, **kwargs):
        super().__init__()
        self.default_pars(
            init_prob=ss.bernoulli(p=0.9),  # Probability that a newly diagnosed person will initiate treatment
            future_coverage={'year': 2022, 'prop': 0.85},
        )
        self.update_pars(pars, **kwargs)
        self.coverage_data = coverage_data
        self.coverage = None  # Set below
        self.coverage_format = None  # Set below
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        data = self.coverage_data
        if data is not None:
            if (len(data.columns) > 1) or (data.columns[0] not in ['n_art', 'p_art']):
                errormsg = 'Expecting a dataframe with a single column labeled n_art or p_art'
                raise ValueError(errormsg)
            colname = data.columns[0]
            self.coverage_format = colname
            self.coverage = sc.smoothinterp(sim.yearvec, data.index.values, data[colname].values)
        self.initialized = True
        return

    def apply(self, sim):
        """
        Apply the ART intervention at each time step. Put agents on and off ART and adjust based on data.
        """
        hiv = sim.diseases.hiv
        inf_uids = hiv.infected.uids

        # Figure out how many people should be treated
        if sim.year < self.pars.future_coverage['year']:
            if self.coverage is None:
                n_to_treat = 0
            else:
                if self.coverage_format == 'n_art':
                    n_to_treat = int(self.coverage[sim.ti]/sim.pars.pop_scale)
                elif self.coverage_format == 'p_art':
                    n_to_treat = int(self.coverage[sim.ti]*len(inf_uids))
        else:
            p_cov = self.pars.future_coverage['prop']
            n_to_treat = int(p_cov*len(inf_uids))

        # Firstly, check who is stopping ART
        if hiv.on_art.any():
            stopping = hiv.on_art & (hiv.ti_stop_art <= sim.ti)
            if stopping.any():
                hiv.stop_art(stopping.uids)

        # Next, see how many people we need to treat vs how many are already being treated
        on_art = hiv.on_art

        # A proportion of newly diagnosed agents onto ART will be willing to initiate ART
        diagnosed = hiv.ti_diagnosed == sim.ti
        if len(diagnosed.uids):
            dx_to_treat = self.pars.init_prob.filter(diagnosed.uids)

            # Figure out if there are treatment spots available and if so, prioritize newly diagnosed agents
            n_available_spots = n_to_treat - len(on_art.uids)
            if n_available_spots > 0:
                self.prioritize_art(sim, n=n_available_spots, awaiting_art_uids=dx_to_treat)

        # Apply correction to match ART coverage data:
        self.art_coverage_correction(sim, target_coverage=n_to_treat)

        # Adjust rel_sus for protected unborn agents
        if hiv.on_art[sim.people.pregnancy.pregnant].any():
            mother_uids = (hiv.on_art & sim.people.pregnancy.pregnant).uids
            infants = sim.networks.maternalnet.find_contacts(mother_uids)
            hiv.rel_sus[ss.uids(infants)] = 0

        return

    def prioritize_art(self, sim, n=None, awaiting_art_uids=None):
        """
        Prioritize ART to n agents among those awaiting treatment
        """
        hiv = sim.diseases.hiv
        if awaiting_art_uids is None:
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids

        # Enough spots for everyone
        if n > len(awaiting_art_uids):
            start_uids = awaiting_art_uids

        # Not enough spots - construct weights based on CD4 count and care seeking
        else:
            cd4_counts = hiv.cd4[awaiting_art_uids]
            care_seeking = hiv.care_seeking[awaiting_art_uids]
            weights = cd4_counts*(1/care_seeking)
            choices = np.argsort(weights)[:n]
            start_uids = awaiting_art_uids[choices]

        hiv.start_art(start_uids)

        return

    def art_coverage_correction(self, sim, target_coverage=None):
        """
        Adjust ART coverage to match data
        """
        hiv = sim.diseases.hiv
        on_art = hiv.on_art

        # Too many agents on treatment -> remove
        if len(on_art.uids) > target_coverage:

            # Agents with the highest CD4 counts will go off ART:
            n_to_stop = int(len(on_art.uids) - target_coverage)
            on_art_uids = on_art.uids

            # Construct weights and choice distribution
            cd4_counts = hiv.cd4[on_art_uids]
            care_seeking = hiv.care_seeking[on_art_uids]
            weights = cd4_counts/care_seeking
            choices = np.argsort(-weights)[:n_to_stop]
            stop_uids = on_art_uids[choices]

            hiv.ti_stop_art[stop_uids] = sim.ti
            hiv.stop_art(stop_uids)

        # Not enough agents on treatment -> add
        elif len(on_art.uids) < target_coverage:
            n_to_add = target_coverage - len(on_art.uids)
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids
            self.prioritize_art(sim, n=n_to_add, awaiting_art_uids=awaiting_art_uids)


class VMMC(ss.Intervention):
    def __init__(self, pars=None, coverage_data=None, eligibility=None, **kwargs):
        super().__init__()
        self.default_pars(
            future_coverage={'year': 2022, 'prop': 0.1},
            eff_circ = 0.6,  # Evidence of a 60% reduction in risk of HIV acquisition: https://www.who.int/teams/global-hiv-hepatitis-and-stis-programmes/hiv/prevention/voluntary-medical-male-circumcision
        )
        self.update_pars(pars, **kwargs)

        # Coverage data - can be number or proportion
        self.coverage_data = coverage_data
        self.coverage = None  # Set below
        self.coverage_format = None  # Set below
        self.eligibility = eligibility  # Determines denominator for coverage if given as a proportion

        # States
        self.willingness = ss.FloatArr('willingness', default=ss.random())  # Willingness to undergo VMMC
        self.circumcised = ss.BoolArr('circumcised', default=False)
        self.ti_circumcised = ss.FloatArr('ti_circumcised')

        return

    def init_pre(self, sim):
        super().init_pre(sim)

        # Handle coverage dataa
        data = self.coverage_data
        if data is not None:
            if (len(data.columns) > 1) or (data.columns[0] not in ['n_vmmc', 'p_vmmc']):
                errormsg = 'Expecting a dataframe with a single column labeled n_vmmc or p_vmmc'
                raise ValueError(errormsg)
            colname = data.columns[0]
            self.coverage_format = colname
            self.coverage = sc.smoothinterp(sim.yearvec, data.index.values, data[colname].values)

        self.init_results()

        return

    def init_post(self):
        super().init_post()
        return

    def init_results(self):
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_circumcisions', npts, dtype=float, scale=True, label="New circumcisions"),
            ss.Result(self.name, 'n_circumcised', npts, dtype=float, scale=True, label="Number circumcised")
        ]
        return

    def apply(self, sim):
        m_uids = sim.people.male.uids

        # Figure out how many people should be circumcised
        if sim.year < self.pars.future_coverage['year']:
            if self.coverage is None:
                n_to_circ = 0
            else:
                if self.coverage_format == 'n_vmmc':
                    n_to_circ = int(sim.dt*self.coverage[sim.ti]/sim.pars.pop_scale)
                elif self.coverage_format == 'p_vmmc':
                    n_to_circ = int(sim.dt*self.coverage[sim.ti]*len(m_uids))
        else:
            p_cov = self.pars.future_coverage['prop']
            n_to_circ = int(sim.dt*p_cov*len(m_uids))

        if n_to_circ > 0:
            # Find who's eligible to circumcise
            eligible_uids = (sim.people.male & ~self.circumcised).uids  # Apply eligiblity
            weights = self.willingness[eligible_uids]
            choices = np.argsort(-weights)[:n_to_circ]
            new_circs = eligible_uids[choices]

            self.circumcised[new_circs] = True
            self.ti_circumcised[new_circs] = sim.ti

        self.results['new_circumcisions'][sim.ti] = n_to_circ
        self.results['n_circumcised'][sim.ti] = count(self.circumcised)

        # Reduce rel_sus
        sim.diseases.hiv.rel_sus[self.circumcised] *= 1-self.pars.eff_circ

        return


# %% Gonorrhea classes
__all__ += ["GonorrheaTreatment", "UpdateDrugs"]


class GonorrheaTreatment(STITreatment):
    """
    Treatment for gonorrhea infection.
        - successful treatment clears infection immediately
        - unsuccessful treatment reduces dur_inf and rel_trans, but results in lower rel_treat
        - unnecessary treatment results in lower rel_treat
    """
    def __init__(self, pars=None, eligibility=None, max_capacity=None, years=None, *args, **kwargs):
        super().__init__(disease='ng', eligibility=eligibility, max_capacity=max_capacity, years=years, *args)
        self.requires = ['ng', 'structuredsexual']
        self.default_pars(
            base_treat_eff=0.96,
            treat_eff=ss.bernoulli(p=0),  # Reset each time step depending on base_treat_eff and population AMR
            rel_treat_unsucc=0.01,
            rel_treat_unneed=0.005,
        )
        self.update_pars(pars, **kwargs)

        # States
        self.add_states(
            ss.FloatArr('rel_treat', default=1),  # How well a person will respond to treatment
        )

    def init_post(self):
        super().init_post()
        results = [
            ss.Result('ng', 'rel_treat', self.sim.npts, dtype=float, scale=False),
        ]
        self.results += results
        self.sim.diseases.ng.results += results
        return

    @staticmethod
    def change_states(sim, treat_succ):
        """ Change the states of people who are treated """
        sim.diseases.ng.clear_infection(treat_succ)
        sim.diseases.ng.wipe_dates(treat_succ)

    def set_treat_eff(self, uids):
        new_treat_eff = self.rel_treat[uids] * self.pars.base_treat_eff
        self.pars.treat_eff.set(new_treat_eff)
        return

    def apply(self, sim):
        """
        Apply treatment. On each timestep, this method will add eligible people who are willing to accept treatment to a
        queue, and then will treat as many people in the queue as there is capacity for.
        """
        treat_uids = super().apply(sim)

        # Change treatment resistance for those unsuccessfully treated
        treat_unsucc = self.outcomes['unsuccessful']
        if len(treat_unsucc):
            self.rel_treat[treat_unsucc] *= (1 - self.pars.rel_treat_unsucc)
        treat_unneed = self.outcomes['unnecessary']
        if len(treat_unneed):
            self.rel_treat[treat_unneed] *= (1 - self.pars.rel_treat_unneed)

        return treat_uids

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        treat_uids = (self.ti_treated == ti).uids
        # Add mean rel_treat among those newly treated
        if len(treat_uids):
            self.sim.diseases.ng.results['rel_treat'][ti] = np.mean(self.rel_treat[treat_uids])

        return


class UpdateDrugs(ss.Intervention):
    """
    An intervention that samples rel_treat and updates the rel_treat values if they fall below
    a given level.
    """
    def __init__(self, pars=None, eligibility=None, years=None, *args, **kwargs):
        super().__init__(*args)
        self.requires = ['ng', 'gonorrheatreatment']
        self.default_pars(
            threshold_amr=0.05
        )
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        self.years = years
        self.add_states(
            ss.BoolArr('rel_treat_prev'),  # Store a copy of AMR to the previous regimen
        )
        self.change_time = None

    def init_post(self):
        super().init_post()
        self.results += [
        ]

    def apply(self, sim):
        target_uids = self.check_eligibility(sim)
        pop_rel_treat = np.mean(self.sim.people.rel_treat[target_uids])
        if pop_rel_treat < self.pars.threshold_amr:
            self.sim.people.rel_treat_prev = sc.dcp(self.sim.people.rel_treat)
            self.sim.people.rel_treat[:] = 1  # Reset
            self.change_time = self.sim.year
        return
