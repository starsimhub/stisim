"""
Define interventions for STIsim
"""

import starsim as ss
import numpy as np
import pandas as pd
from collections import defaultdict
import sciris as sc
import functools
import itertools

__all__ = ['HIVTest', 'ART', 'VMMC', "PartnerNotification", "SyphVaccine"]


class HIVTest(ss.Intervention):
    """
    Base class for HIV testing

    Args:
         prob           (float/arr)     : annual probability of eligible people being tested
         eligibility    (inds/callable) : indices OR callable that returns inds
         label          (str)           : the name of screening strategy
         kwargs         (dict)          : passed to Intervention()
    """

    def __init__(self, pars=None, test_prob_data=None, years=None, start_year=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.default_pars(
            rel_test=1,
        )
        self.update_pars(pars, **kwargs)

        # Set testing probabilities and years
        if years is not None and start_year is not None:
            errormsg = 'Provide either years or start_year, not both.'
            raise ValueError(errormsg)
        self.years = years
        self.start_year = start_year
        self.test_prob_data = test_prob_data
        self.test_prob = ss.bernoulli(self.make_test_prob_fn)

        # Set eligibility
        self.eligibility = eligibility
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed

        # States
        self.tested = ss.BoolArr('tested', default=False)
        self.ti_tested = ss.FloatArr('ti_tested')
        self.diagnosed = ss.BoolArr('diagnosed', default=False)
        self.ti_diagnosed = ss.FloatArr('ti_diagnosed')

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start_year is None:
            self.start_year = self.years[0]
        self.init_results()
        return

    def init_results(self):
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_diagnoses', npts, dtype=float, scale=True),
            ss.Result(self.name, 'new_tests', npts, dtype=int, scale=True)]
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

    def apply(self, sim):
        if sim.year > self.start_year:
            hiv = sim.diseases.hiv

            # Find who's eligible to test, who gets a test, and who is diagnosed
            eligible_uids = self.check_eligibility(sim)  # Apply eligiblity
            if len(eligible_uids):
                tester_uids = self.test_prob.filter(eligible_uids)
                if len(tester_uids):
                    # Add results and states for testers
                    self.results['new_tests'][sim.ti] += len(tester_uids)
                    self.tested[tester_uids] = True
                    self.ti_tested[tester_uids] = sim.ti

                    # Add results and states for diagnoses
                    pos_uids = tester_uids[hiv.infected[tester_uids]]
                    self.results['new_diagnoses'][sim.ti] += len(pos_uids)
                    self.diagnosed[pos_uids] = True
                    self.ti_diagnosed[pos_uids] = sim.ti
                    hiv.diagnosed[pos_uids] = True
                    hiv.ti_diagnosed[pos_uids] = sim.ti

        return

    def check_eligibility(self, sim):
        if self.eligibility is not None:
            uids = self.eligibility(sim).uids
        else:
            uids = sim.people.alive.uids
        return uids


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
            ss.Result(self.name, 'new_circumcisions', npts, dtype=float, scale=True),
            ss.Result(self.name, 'n_circumcised', npts, dtype=float, scale=True)]
        return

    def apply(self, sim):
        hiv = sim.diseases.hiv
        males = sim.people.male
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
        self.results['n_circumcised'][sim.ti] = np.count_nonzero(self.circumcised)

        # Reduce rel_sus
        sim.diseases.hiv.rel_sus[self.circumcised] *= 1-self.pars.eff_circ

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


class SyphVaccine(ss.Intervention):
    def __init__(self, pars=None, years=None, start_year=None, eligibility=None, target_coverage=None, name=None, label=None, 
                  **kwargs):
        super().__init__(name=name, label=label)
        self.name = 'syph_vaccine'
        self.requires = 'syphilis'
        self.default_pars(
            efficacy=0.9,
            dur_protection=12, # half-life of exponential decay
            dur_reach_peak=2, # 2 months until efficacy is reached
            immunity_init=ss.uniform(low=0.7, high=0.9),
            nab_boost_infection=0.9, # Multiply base immunity by this factor. 1=no change in immunity, 0=full immunity, no reinfection
            nab_boost_vaccination=0.7, # Multiply base immunity by this factor. 1=no change in immunity, 0=full immunity, no reinfection
            prevent_infection=0.05,
            prevent_transmission_susceptible=0,
            prevent_transmission_primary=0.9,
            prevent_transmission_secondary=0.1,
            prevent_transmission_tertiary=0.05,
            prevent_transmission_latent=0.05,
            reduce_dur_primary=0.7,
            reduce_dur_secondary=0.2,
            reduce_p_reactivate=0.5
        )
        self.update_pars(pars, **kwargs)
        
        # Years
        if years is not None and start_year is not None:
            errormsg = 'Provide either years or start_year, not both.'
            raise ValueError(errormsg)
        self.years = years
        self.start_year = start_year

        # Eligibility
        self.eligibility = eligibility

        # States
        self.add_states(
            ss.BoolArr('vaccinated'),
            ss.FloatArr('ti_vaccinated'),
            ss.FloatArr('doses', default=0),
            ss.FloatArr('immunity_trans', default=1),
            ss.FloatArr('immunity_inf', default=1),
            ss.FloatArr('ti_nab_event'),
        )

        # Vaccine 
        self.current_coverage = 0
        self.target_coverage = 0.7
        self.dose_interval = None
        self._immunity_timecourse = None
        self._protection_timecourse = None

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start_year is None:
            self.start_year = self.years[0]
        self.init_results()
        self._immunity_timecourse = self.immunity_timecourse()
        self._protection_timecourse = self.immunity_timecourse() # For now, assume protection timecourse and immunity timecourse are equal 
        return

    def init_results(self):
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_vaccinations', npts, dtype=float, scale=True),
            ss.Result(self.name, 'n_vaccinated', npts, dtype=int, scale=True),
        ]
        return

    def linear_increase(self, length, init_val, slope):
        ''' 
        Calculate linear decay 
        '''
        result = slope * np.ones(length)
        result[0] = init_val
        return result.cumsum()

    def exp_decay(self, t, init_val, half_life):
        '''
        Returns an array of length t with values for the immunity at each time step
        '''
        decay_rate = np.log(2) / half_life if ~np.isnan(half_life) else 0.
        result = init_val * np.exp(-decay_rate * t, dtype=ss.dtypes.float)
        return result
    
    def immunity_timecourse(self):
        
        # Efficacy will increase linearly to its peak value
        linear_increase = self.linear_increase(length=self.pars.dur_reach_peak, init_val=0, slope=self.pars.efficacy/self.pars.dur_reach_peak)
        # Efficacy will then drop exponentially, with half-time corresponding to the duration of proection
        exp_decay = self.exp_decay(t=np.arange(0, self.sim.npts - len(linear_increase)), init_val=self.pars.efficacy, half_life=self.pars.dur_protection)
        # Combine to one array
        timecourse = np.concatenate([linear_increase, exp_decay])
        
        # Return derivative 
        timecourse_derivative = np.diff(timecourse)
        return timecourse_derivative

    def check_eligibility(self, sim):
        if self.eligibility is not None:
            uids = self.eligibility(sim).uids
        else:
            uids = sim.people.alive.uids
        return uids

    def get_targets(self, sim):
        target_uids = ss.uids()
        eligible_uids = self.check_eligibility(sim)  # Apply eligiblity
        current_vaccinated = self.vaccinated.uids

        n_current_vaccinated = len(current_vaccinated)
        n_target_vaccinated = len(eligible_uids) * self.target_coverage
        n_to_vaccinate = int(n_target_vaccinated - n_current_vaccinated)
       
        if n_to_vaccinate > 0:
            # Pick eligible, non-vaccinated agents randomly
            # eligible_uids = eligible_uids & (~self.vaccinated).uids # Allow for multiple doses
            bools = ss.random(strict=False).rvs(len(eligible_uids))
            choices = np.argsort(bools)[:n_to_vaccinate]
            target_uids = eligible_uids[choices]

        return target_uids

    def update_natural_immunity(self, sim):
        """
        Update base immunity for individuals, who got infected at this timestep
        """
        # Extract parameters and indices
        syph = sim.diseases.syphilis
        new_syphilis = syph.ti_infected == sim.ti
        has_nabs = self.ti_nab_event.notnan

        prior_nab_uids = (new_syphilis & has_nabs).uids
        no_prior_nab_uids = (new_syphilis & ~has_nabs).uids

        # 1) Individuals that already have NAbs from a previous vaccination/infection have their Immunity levels boosted
        # immunity_inf = 1 -> No Immunity, immunity_inf = 0 -> full immunity, no reinfection
        if len(prior_nab_uids):
            boost_factor = self.pars.nab_boost_infection
            self.immunity_inf[prior_nab_uids] *= boost_factor

        # 2) Individuals without prior NAbs are assigned an initial level drawn from a distribution.
        if len(no_prior_nab_uids):
            self.immunity_inf[no_prior_nab_uids] = self.pars.immunity_init.rvs(no_prior_nab_uids)

        # Ensure values between 0 and 1
        # self.immunity_inf = np.where(self.immunity_inf >= 0, self.immunity_inf, 0)  # Make sure immunity doesn't drop below 0
        # self.immunity_inf = np.where(self.immunity_inf <= 1, self.immunity_inf, 1)  # Make sure immunity doesn't exceed 1

        # Update time of NAb event
        self.ti_nab_event[new_syphilis] = sim.ti
        return

    def vaccinate(self, sim, uids, update_immunity_by_vaccination=True):
        """
        Vaccinate
        """
        # Set states
        self.vaccinated[uids] = True
        self.ti_vaccinated[uids] = sim.ti
        self.ti_nab_event[uids] = sim.ti
        # Update number of doses
        self.doses[uids] += 1
        if update_immunity_by_vaccination:
            # Boost Immunity for >1 dose
            boost_uids = uids[self.doses[uids] > 1]
            if len(boost_uids):
                self.boost_immunity_by_vaccination(sim, boost_uids)
            # Update Immunity
            self.update_immunity_by_vaccination(sim)

    def boost_immunity_by_vaccination(self, sim, uids):
        """
        Boost immunity level for newly vaccinated individuals, who received their second, third, etc. dose
        """
        boost_factor = self.pars.nab_boost_vaccination
        self.immunity_inf[uids] *= boost_factor
        return
        
    def update_immunity_by_vaccination(self, sim):
        """
        Update Immunity levels for both vaccinated and unvaccinated individuals 
        """
        syph = sim.diseases.syphilis

        # Vaccinated Individuals
        vaccinated_bool = self.vaccinated
        vaccinated_uids = vaccinated_bool.uids
        
        # Update protection against infection and transmission for vaccinated inidivuals 
        # Do this by infection state to differentiate between different transmission reductions per state  
        for state in ['susceptible', 'primary', 'secondary', 'tertiary', 'latent']:
            state = f'{state}'
            state_bools = getattr(sim.diseases.syphilis, state)
            uids = vaccinated_uids & state_bools.uids
            if len(uids):
                ti_since_boost_vaccinated = sim.ti - self.ti_nab_event[uids].astype(ss.dtypes.int)
                ti_since_boost_vaccinated = np.minimum(ti_since_boost_vaccinated, len(self._protection_timecourse) - 1).astype(ss.dtypes.int)  # Max out protection
                
                immunity = self._immunity_timecourse[ti_since_boost_vaccinated]
                protection = self._protection_timecourse[ti_since_boost_vaccinated]
                prevent_transmission_param = self.pars[f'prevent_transmission_{state}']
                
                # Update protection against infection for vaccinated agents
                self.immunity_inf[uids] -= self.pars.prevent_infection * immunity

                # Update protection against transmission for vaccinated inidivuals
                self.immunity_trans[uids] -= prevent_transmission_param * protection

        # Ensure values between 0 and 1
        # self.immunity_inf = np.where(self.immunity_inf >= 0, self.immunity_inf, 0)  # Make sure immunity doesn't drop below 0
        # self.immunity_inf = np.where(self.immunity_inf <= 1, self.immunity_inf, 1)  # Make sure immunity doesn't exceed 1

        # Set rel trans and rel sus
        rel_trans, rel_sus = self.compute_trans_sus(sim)
        syph.rel_trans.set(uids=sim.people.auids, new_vals = rel_trans)
        syph.rel_sus.set(uids=sim.people.auids, new_vals = rel_sus)

        # TODO Update dur_primary, dur_secondary, p_reinfection
        # ti_primary
        # ti_secondary
        return

    def compute_trans_sus(self, sim):
        syph = sim.diseases.syphilis
        rel_trans = syph.rel_trans * syph.infectious * self.immunity_trans
        rel_sus = syph.rel_sus * syph.susceptible * self.immunity_inf
        return rel_trans, rel_sus

    def apply(self, sim):
        syph = sim.diseases.syphilis
        self.update_natural_immunity(sim)
        if sim.year > self.start_year:
            target_uids = self.get_targets(sim)
            if len(target_uids):
                self.vaccinate(sim, target_uids)
                self.results['new_vaccinations'][sim.ti] += len(target_uids)
            else:
                self.update_immunity_by_vaccination(sim)

        return
