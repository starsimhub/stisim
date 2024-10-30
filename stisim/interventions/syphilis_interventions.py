"""
Define syphilis interventions for STIsim
"""

import starsim as ss
import numpy as np
import pandas as pd
import sciris as sc
from stisim.interventions.base_interventions import STIDx, STITest, STITreatment
from stisim.utils import TimeSeries
from sciris import randround as rr  # Since used frequently

__all__ = ["SyphDx", "SyphTx", "NewbornTreatment", "SyphTest", "ANCSyphTest", "NewbornSyphTest", "SyphVaccine"]


class SyphDx(STIDx):
    def __init__(self, df, *args, **kwargs):
        super().__init__(df, 'syphilis', *args, **kwargs)
        return


class SyphTx(STITreatment):
    """
    Treat a fixed number of people each timestep.
    """

    def __init__(self, pars=None, max_capacity=None, years=None, eligibility=None, name=None, **kwargs):
        super().__init__(diseases='syphilis', name=name, eligibility=eligibility, years=years, max_capacity=max_capacity)
        self.define_pars(
            rel_treat_prob=1,
            treat_prob=ss.bernoulli(p=1),
            treat_eff=ss.bernoulli(p=0.95),
            fetus_age_cutoff_treat_eff=-0.25,  # Reduced treatment efficacy for fetuses in the last trimester
            treat_eff_reduced=ss.bernoulli(p=0.2)  # Reduced efficacy for fetuses older than cut off
        )
        self.update_pars(pars, **kwargs)
        return

    def administer(self, sim, uids, disease='syphilis', return_format='dict'):
        """ Administer treatment, keeping track of unnecessarily treated individuals """

        inf = sim.diseases[disease].infected
        sus = sim.diseases[disease].susceptible
        inf_uids = uids[inf[uids]]
        sus_uids = uids[sus[uids]]

        successful = self.pars.treat_eff.filter(inf_uids)
        unsuccessful = np.setdiff1d(inf_uids, successful)
        unnecessary = sus_uids

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
        elif return_format == 'array':
            output = successful

        return output

    def change_states(self, disease, treat_succ):
        """ Change the states of people who are treated """
        self.sim.diseases[disease].primary[treat_succ] = False
        self.sim.diseases[disease].secondary[treat_succ] = False
        self.sim.diseases[disease].latent[treat_succ] = False
        self.sim.diseases[disease].tertiary[treat_succ] = False
        self.sim.diseases[disease].ti_primary[treat_succ] = np.nan
        self.sim.diseases[disease].ti_secondary[treat_succ] = np.nan
        self.sim.diseases[disease].ti_latent[treat_succ] = np.nan
        self.sim.diseases[disease].ti_tertiary[treat_succ] = np.nan
        self.sim.diseases[disease].susceptible[treat_succ] = True
        self.sim.diseases[disease].infected[treat_succ] = False

    def treat_fetus(self, sim, mother_uids):
        """
        Treat fetuses of successfully treated mothers.
        Birth outcomes of successfully treated fetuses will be updated. If fetus is not yet infected, the treatment doesn't do anything.
        Treatment success depends on when in the pregnancy the mother gets treated
        with reduced treatment efficacy above the cut off age
        """
        # Get Fetuses
        children_uids = ss.uids(sim.networks.maternalnet.find_contacts(mother_uids))  # Children and Fetuses of treated mothers
        fetus_gt_cutoff = children_uids[sim.people.age[children_uids] > self.pars.fetus_age_cutoff_treat_eff]  # Fetuses above cutoff age
        fetus_st_cutoff = children_uids[sim.people.age[children_uids] <= self.pars.fetus_age_cutoff_treat_eff]  # Fetuses below cutoff age

        # Treat fetuses above cutoff age (e.g. in the last trimester) with a reduced treatment efficacy
        treat_gt_cutoff = self.pars.treat_eff_reduced.filter(fetus_gt_cutoff)

        # Treat fetuses below cutoff age (e.g. in the first two trimesters)
        treat_st_cutoff = self.pars.treat_eff.filter(fetus_st_cutoff)

        # Combine
        treat_uids = treat_gt_cutoff.concat(treat_st_cutoff)

        # Store results - Failure
        sim.diseases.syphilis.results['new_fetus_treated_failure'][sim.ti] = len(fetus_gt_cutoff) + len(fetus_st_cutoff) - len(treat_uids)
        # Store results - Unnecessary. Subtract successful uids below
        sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][sim.ti] += len(treat_uids)

        # Change birth outcomes for successfully treated unborn babies
        # If fetus is not yet infected, this will do nothing
        if len(treat_uids):
            # Birth outcomes must be modified to add probability of susceptible birth
            birth_outcomes = sim.diseases.syphilis.pars.birth_outcome_keys

            # Change states
            for oi, outcome in enumerate(birth_outcomes):
                ti_outcome = f'ti_{outcome}'
                vals = getattr(sim.diseases.syphilis, ti_outcome)
                successful_uids = treat_uids & vals.notnan.uids
                vals[successful_uids] = np.nan
                setattr(sim.diseases.syphilis, ti_outcome, vals)

                # Store results - Success
                sim.diseases.syphilis.results['new_fetus_treated_success'][sim.ti] += len(successful_uids)
                # Store results - Unnecessary. Subtract the uids that were successful
                sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][sim.ti] -= len(successful_uids)

        return

    def step(self):
        """
        Apply treatment
        """
        sim = self.sim
        treat_uids = super().step()
        # Treat unborn babies of successfully treated mothers
        treat_pregnant_uids = sim.people.pregnancy.pregnant.uids & self.outcomes['successful']
        if len(treat_pregnant_uids):
            self.treat_fetus(sim, mother_uids=treat_pregnant_uids)
        return treat_uids


class NewbornTreatment(SyphTx):

    def change_states(self, disease, treat_succ):
        """ Change states of congenital cases """
        self.sim.diseases[disease].congenital[treat_succ] = False
        self.sim.diseases[disease].ti_congenital[treat_succ] = np.nan
        # sim.diseases.syphilis.susceptible[treat_succ] = True  # Leave this out for now

    def administer(self, sim, uids, disease, return_format='dict'):
        """ Administer treatment to newborns """
        sus = sim.diseases[disease].susceptible
        con = sim.diseases[disease].congenital
        sus_uids = uids[sus[uids]]
        con_uids = uids[con[uids]]

        successful = self.pars.treat_eff.filter(con_uids)
        unsuccessful = np.setdiff1d(con_uids, successful)
        unnecessary = sus_uids

        # Return outcomes
        if return_format == 'dict':
            output = {'successful': successful, 'unsuccessful': unsuccessful, 'unnecessary': unnecessary}
        elif return_format == 'array':
            output = successful

        return output


class SyphTest(STITest):
    """ Base class for syphilis tests """

    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None, name=None, label=None, newborn_test=None,
                 **kwargs):
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, eligibility=eligibility, product=product, name=name, label=label, **kwargs)
        self.define_pars(
            linked=True,
        )
        self.update_pars(pars, **kwargs)
        # Store optional newborn test intervention
        self.newborn_test = newborn_test
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.test_prob_data = self.process_data(sim)
        return

    def process_data(self, sim):
        """ Turn dataframe into a dictionary """
        if isinstance(self.test_prob_data, pd.DataFrame):
            df = self.test_prob_data.set_index(['risk_group', 'sex', 'sw'])
            df = df.pivot(columns='year', values='symp_test_prob')
            dd = df.to_dict(orient='index')
            for group, vals in dd.items():
                dd[group] = sc.smoothinterp(sim.timevec, list(vals.keys()), list(vals.values()), smoothness=0)
            return dd
        else:
            return self.test_prob_data

    @staticmethod
    def make_test_prob_fn(self, sim, uids):
        """ Process symptomatic testing probabilites over time by sex and risk group """

        if sc.isnumber(self.test_prob_data):
            test_prob = self.test_prob_data
        elif isinstance(self.test_prob_data, TimeSeries):
            test_prob = self.test_prob_data.interpolate(sim.now)
        elif sc.checktype(self.test_prob_data, 'arraylike'):
            year_ind = sc.findnearest(self.years, sim.now)
            test_prob = self.test_prob_data[year_ind]
        elif isinstance(self.test_prob_data, dict):
            test_prob = pd.Series(index=uids)
            n_risk_groups = sim.networks.structuredsexual.pars.n_risk_groups
            for rg in range(n_risk_groups):
                for sex in ['female', 'male']:
                    for sw in [0, 1]:
                        dkey = (rg, sex, sw)
                        conditions = (sim.people[sex] & (sim.networks.structuredsexual.risk_group == rg))
                        if sw:
                            if sex == 'female': conditions = conditions & sim.networks.structuredsexual.fsw
                            if sex == 'male':   conditions = conditions & sim.networks.structuredsexual.client
                        test_prob[conditions[uids]] = self.test_prob_data[dkey][sim.ti]
        else:
            errormsg = 'Format of test_prob_data must be float, array, or dict.'
            raise ValueError(errormsg)

        # Scale and validate
        test_prob = test_prob * self.pars.rel_test
        if not self.pars.linked:
            test_prob = test_prob * self.dt
        test_prob = np.clip(test_prob, a_min=0, a_max=1)

        return test_prob

    def step(self, uids=None):
        super().step(uids=uids)
        sim = self.sim
        if (sim.now >= self.start) & (sim.now < self.stop):
            # Schedule newborn tests if the mother is positive
            if self.newborn_test is not None:
                new_pos = self.ti_positive == self.sim.ti
                if new_pos.any():
                    pos_mother_inds = np.in1d(sim.networks.maternalnet.p1, new_pos.uids)
                    unborn_uids = sim.networks.maternalnet.p2[pos_mother_inds]
                    ti_births = sim.networks.maternalnet.edges.end[pos_mother_inds].astype(int)
                    self.newborn_test.schedule(unborn_uids, ti_births)

        return

    def update_results(self):
        ti = self.sim.ti
        super().update_results()
        if 'positive' in self.outcomes.keys() and (self.ti_positive == ti).any():
            self.update_positives()
        if 'negative' in self.outcomes.keys() and (self.ti_negative == ti).any():
            self.update_negatives()

        return

    def update_positives(self):
        ti = self.sim.ti
        new_pos = self.ti_positive == ti
        # Count true/false positives
        false_pos = np.count_nonzero(self.sim.diseases.syphilis.susceptible[new_pos])
        true_pos = np.count_nonzero(self.sim.diseases.syphilis.infected[new_pos])
        self.sim.diseases.syphilis.results['new_false_pos'][ti] += false_pos
        self.sim.diseases.syphilis.results['new_true_pos'][ti] += true_pos
        return

    def update_negatives(self):
        ti = self.sim.ti
        new_neg = self.ti_negative == ti
        # Count true/false negatives
        false_neg = np.count_nonzero(self.sim.diseases.syphilis.infected[new_neg])
        true_neg = np.count_nonzero(self.sim.diseases.syphilis.susceptible[new_neg])
        self.sim.diseases.syphilis.results['new_false_neg'][ti] += false_neg
        self.sim.diseases.syphilis.results['new_true_neg'][ti] += true_neg
        return


class ANCSyphTest(SyphTest):
    """
    Test given to pregnant women
    Need to adjust timing using Trivedi (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7138526/)
    """

    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None, name=None, label=None, newborn_test=None,
                 **kwargs):
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, eligibility=eligibility, product=product, name=name, label=label, **kwargs)
        self.define_pars(
            linked=True,
        )
        self.update_pars(pars, **kwargs)
        self.test_timing = ss.randint(1, 9)
        if self.eligibility is None:
            self.eligibility = lambda sim: sim.demographics.pregnancy.pregnant
        return

    def get_testers(self, sim):
        # For ANC testing, only administer scheduled tests
        return (self.ti_scheduled == sim.ti).uids

    def schedule_tests(self):
        """ Schedule a test for newly pregnant women """
        sim = self.sim
        newly_preg = (sim.demographics.pregnancy.ti_pregnant == sim.ti).uids
        self.test_prob.pars['p'] = self.make_test_prob_fn(self, sim, newly_preg)
        will_test = self.test_prob.filter(newly_preg)
        ti_test = sim.ti + self.test_timing.rvs(will_test)
        self.ti_scheduled[will_test] = ti_test

    def step(self):
        self.schedule_tests()  # Check for newly pregnant women so they can be added to the schedule
        return super().step()


class NewbornSyphTest(SyphTest):
    """
    Test given to newborns if the mother was confirmed to have syphilis at any stage of the pregnancy
    """

    def get_testers(self, sim):
        # For newborn testing, only administer scheduled tests, and account for probability of testing at this point
        eligible_uids = (self.ti_scheduled == sim.ti).uids
        accept_uids = self.test_prob.filter(eligible_uids)
        return accept_uids


class SyphVaccine(ss.Intervention):
    def __init__(self, pars=None, years=None, coverage_data=None, start_year=None, eligibility=None, name=None, label=None,
                 **kwargs):
        super().__init__(name=name, label=label)
        self.name = 'syph_vaccine'
        self.requires = 'syphilis'
        self.define_pars(
            # Dose parameters
            daily_num_doses=None,  # Daily supply of vaccine doses (unscaled), None equals unlimited supply
            dose_interval=ss.lognorm_ex(mean=3, std=1 / 12),  # Assume every 3 years for now
            p_second_dose=ss.bernoulli(p=0.6),  # Probability that a a person, who received 1 dose, comes back for a second dose
            p_third_dose=ss.bernoulli(p=0.6),  # Probability that a person, who receieved 2 doses, comes back for a third dose. More likely?

            # Immunity parameters
            # - Efficacy
            first_dose_efficacy=0.5,  # vaccine efficacy applied to general pop
            second_dose_efficacy=0.6,
            third_dose_efficacy=0.75,
            rel_efficacy_red_maternal=0.2,  # Relative reduction in efficacy for maternal network, 20% reduction of 90% efficacy -> efficacy of 72%
            # - Protection
            dur_protection=5,  # in years, expected duration of protection, half-life of exponential decay
            dur_reach_peak=0.5,  # in years, Assume 6 months until efficacy is reached
            # - Update immunity
            immunity_init=ss.lognorm_ex(mean=0.001, std=0.001),
            nab_boost_infection=0.05,  # Multiply base immunity by this factor. 1=no change in immunity, 0=full immunity, no reinfection

            prevent_infection=0.05,
            prevent_transmission_susceptible=0.95,
            prevent_transmission_primary=0.95,
            prevent_transmission_secondary=0.95,
            prevent_transmission_tertiary=0.05,
            prevent_transmission_latent=0.05,
            # - Update duration of infection
            reduce_dur_primary=0.2,
            reduce_dur_secondary=0.2,
            # - Reduce probability to reactive secondary state
            reduce_p_reactivate=0.5  # Reduce p_reactivate in syphilis disease module by this factor
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
        self.define_states(
            ss.BoolArr('vaccinated'),
            ss.FloatArr('immunity', default=0),
            ss.FloatArr('rel_trans_immunity', default=0),
            ss.FloatArr('rel_trans_immunity_maternal', default=0),
            ss.FloatArr('rel_sus_immunity', default=0),
            ss.FloatArr('peak_immunity', default=0),
            ss.FloatArr('linear_boost'),
            ss.BoolArr('dur_inf_updated'),
            ss.FloatArr('doses', default=0),
            ss.FloatArr('ti_vaccinated'),
            ss.FloatArr('ti_nab_event'),
            ss.FloatArr('ti_second_dose'),
            ss.FloatArr('ti_third_dose'),
            ss.FloatArr('ti_start_waning'),
        )

        # Vaccine
        self.current_coverage = 0
        self.coverage_data = coverage_data
        self.target_coverage = None  # Set below
        self.coverage_format = None  # Set below

        self.dose_interval = None

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start_year is None:
            self.start_year = self.years[0]

        # Initialize target coverage data
        data = self.coverage_data
        if data is not None:
            colname = data.columns[0]
            self.coverage_format = colname
            self.target_coverage = sc.smoothinterp(sim.yearvec, data.index.values, data[colname].values)
        else:
            self.target_coverage = np.ones(sim.npts)

        # Scale number of doses
        if self.pars.daily_num_doses is not None:
            self.num_doses = rr((self.pars.daily_num_doses * 365 * sim.dt_year) / sim.pars.pop_scale)
        else:
            self.num_doses = None

        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result('new_vaccinated', dtype=float, scale=True, label='Newly vaccinated'),
            ss.Result('n_vaccinated', dtype=int, scale=True, label='Number of vaccinated individuals'),
            ss.Result('n_vaccinated_twice', dtype=int, scale=True, label='Number of vaccinated individuals with two doses'),
            ss.Result('n_vaccinated_triple', dtype=int, scale=True, label='Number of vaccinated individuals with three dose'),
            ss.Result('n_doses', dtype=float, scale=True, label='Number of doses administered'),
        ]
        self.define_results(*results)

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

    def logistic_decay(self, t, init_val):
        """
        Return an array of length t with values for the immunity at each time step
        """
        dt = self.sim.dt_year
        decay_rate = 3 * dt  # Decay rate reprsenting how quickly y decays
        D = 0  # Value of immunity when t -> inf
        y = init_val + (D - init_val) / (1 + np.exp(-decay_rate * (t - rr(self.pars.dur_protection / dt))))
        return y

    def get_immunity_timecourse(self, efficacy, dur_reach_peak, dur_protection):
        """
        Get the derivative of the immunity timecourse.
        Immunity will increase linearly to the vaccine's efficacy level and then decrease exponentially.

            Args:
                efficacy: Vaccine Efficacy
                dur_reach_peak: Parameter to describe how long it takes to reach efficacy
                dur_protection: Parameter to describe how long protection lasts. This will be the half-life of the exponential decay

        """
        dt = self.sim.dt_year
        # Efficacy will increase linearly to its peak value
        linear_increase = self.linear_increase(length=rr(dur_reach_peak / dt), init_val=0, slope=efficacy / rr(dur_reach_peak / dt))
        # Efficacy will then drop exponentially, with half-time corresponding to the duration of protection
        exp_decay = self.exp_decay(t=np.arange(0, self.sim.npts - len(linear_increase)), init_val=efficacy, half_life=rr(dur_protection / dt))
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

    def prioritize_agents(self, sim, num_doses):
        """
        Pick agents to receive a second and third dose. Currently, this is done randomly up to the number of available doses.
        Other options would be to prioritize agents by wait time, or by number of dose (e.g. administer second doses first, then
        third doses, etc.)
        """
        # Eligible for a second and third dose
        eligible_second_dose = self.ti_second_dose <= sim.ti
        eligible_third_dose = self.ti_third_dose <= sim.ti

        # Pick agents randomly
        eligible_uids = eligible_second_dose.uids & eligible_third_dose.uids
        bools = ss.random(strict=False).rvs(len(eligible_uids))
        choices = np.argsort(bools)[:num_doses]
        uids_to_revaccinate = eligible_uids[choices]

        # Pick agents by wait time
        # wait_times_second_dose = sim.ti - self.ti_second_dose[eligible_second_dose.uids]
        # wait_times_third_dose = sim.ti - self.ti_third_dose[eligible_third_dose.uids]
        # wait_times_combined = np.concatenate([wait_times_second_dose, wait_times_third_dose])
        # uids_combined = eligible_second_dose.uids.concat(eligible_third_dose.uids)
        #
        # choices = np.argsort(wait_times_combined)[:remaining_doses]
        # uids_to_revaccinate = uids_combined[choices]

        return uids_to_revaccinate

    def get_targets(self, sim, num_doses=None, target_coverage=None):
        """
        Get uids of agents to get vaccinated this time step.

        1) First, try to reach the target coverage by vaccinating as many eligible agents as possible (up to num_doses).

        2) If target coverage has been reached, and there are doses leftover, administer second and third doses.

        Args:
            num_doses:      available doses at this time step, if None assume unlimited supply
        """
        target_uids = ss.uids()
        eligible_uids = self.check_eligibility(sim)  # Apply eligiblity

        # 1) Use doses to reach target coverage
        current_vaccinated = self.vaccinated.uids
        n_current_vaccinated = len(current_vaccinated)
        n_target_vaccinated = len(eligible_uids) * target_coverage
        # If there unlimited doses (i.e. num_doses=None), vaccinate all, otherwise vaccinate as many as possible
        if num_doses is not None:
            n_to_vaccinate = np.minimum(num_doses, int(n_target_vaccinated - n_current_vaccinated))
        else:
            n_to_vaccinate = int(n_target_vaccinated - n_current_vaccinated)

        target_coverage_uids = ss.uids()
        if n_to_vaccinate > 0:
            # Pick eligible, non-vaccinated agents randomly to each target coverage
            eligible_uids = eligible_uids & (~self.vaccinated).uids
            bools = ss.random(strict=False).rvs(len(eligible_uids))
            choices = np.argsort(bools)[:n_to_vaccinate]
            target_coverage_uids = eligible_uids[choices]

        eligible_second_dose = self.ti_second_dose <= sim.ti
        eligible_third_dose = self.ti_third_dose <= sim.ti

        # 2) Disperse unused doses to agents, who are eligible to receive a second or third dose
        # If there are any unused doses, offer a second and third dose to any vaccinated agents, scheduled to come back for second or third dose
        if num_doses is not None:
            remaining_doses = num_doses - len(target_coverage_uids)
            uids_to_revaccinate = self.prioritize_agents(sim, remaining_doses)
        else:
            uids_to_revaccinate = eligible_second_dose.uids | eligible_third_dose.uids

        # Reset ti_second_dose, and ti_third_dose for agents who received their second and third dose
        get_second_dose_uids = uids_to_revaccinate & eligible_second_dose.uids
        get_third_dose_uids = uids_to_revaccinate & eligible_third_dose.uids
        self.ti_second_dose[get_second_dose_uids] = np.nan
        self.ti_third_dose[get_third_dose_uids] = np.nan

        # Combine all agents that will get vaccinated at this timestep
        target_uids = target_coverage_uids | get_second_dose_uids | get_third_dose_uids

        return target_uids

    def update_natural_immunity(self, sim):
        """
        Update rel_sus_immunity for individuals, who got infected with syphilis at this timestep
        No effect on rel_trans_immunity
        """
        # Extract parameters and indices
        syph = sim.diseases.syphilis
        new_syphilis = syph.ti_infected == sim.ti
        has_nabs = self.ti_nab_event.notnan

        prior_nab_uids = (new_syphilis & has_nabs).uids
        no_prior_nab_uids = (new_syphilis & ~has_nabs).uids

        # 1) Individuals that already have NAbs from a previous vaccination/infection have their Immunity levels boosted
        if len(prior_nab_uids):
            linear_boost = self.pars.nab_boost_infection
            self.immunity[prior_nab_uids] += linear_boost
            self.peak_immunity[prior_nab_uids] = self.immunity[prior_nab_uids]

        # 2) Individuals without prior NAbs are assigned an initial level drawn from a distribution.
        if len(no_prior_nab_uids):
            self.immunity[no_prior_nab_uids] = self.pars.immunity_init.rvs(no_prior_nab_uids)

        # Update time of NAb event
        self.ti_nab_event[new_syphilis] = sim.ti

        # Update peak immunity
        self.peak_immunity[new_syphilis] = np.maximum(self.peak_immunity[new_syphilis], self.immunity[new_syphilis])

        # Start waning in the next time step
        self.ti_start_waning[new_syphilis] = sim.ti + 1

        # Update rel sus and rel trans for new syphilis cases
        self.rel_sus_immunity[new_syphilis] = self.immunity[new_syphilis]
        self.rel_trans_immunity[new_syphilis] = self.immunity[new_syphilis]
        self.rel_trans_immunity_maternal[new_syphilis] = self.immunity[new_syphilis]

        # Wane, only wane unvaccinated agents here
        uids_to_wane = (~self.vaccinated).uids & (self.ti_start_waning <= sim.ti).uids
        self.immunity[uids_to_wane] = self.logistic_decay(sim.ti - self.ti_start_waning[uids_to_wane], self.peak_immunity[uids_to_wane])

        # Update rel sus and rel trans for waning agents
        self.rel_sus_immunity[uids_to_wane] = self.immunity[uids_to_wane]
        self.rel_trans_immunity[uids_to_wane] = self.immunity[uids_to_wane]
        self.rel_trans_immunity_maternal[uids_to_wane] = self.immunity[uids_to_wane]

        return

    def set_linear_boost(self, sim, uids):
        """
        Get linear boost values
        """
        dt = sim.dt_year
        hiv = sim.diseases.hiv
        hiv_pos = hiv.infected.uids
        hiv_pos_uids = uids & hiv_pos
        hiv_neg_uids = uids.remove(hiv_pos_uids)

        # Boost depending on CD4 counts
        cd4_bins = np.array([1000, 500, 350, 200, 50, 0])  # TODO make input
        hiv_rel_linear_boost = np.array([1, 1, 0.9, 0.8, 0.7, 0.5])  # Percentage of HIV negative linear boost

        new_vaccinated_uids = uids[self.ti_vaccinated[uids] == sim.ti]

        # Set this once for newly vaccinated agents
        self.linear_boost[new_vaccinated_uids] = (self.peak_immunity[new_vaccinated_uids] - self.immunity[new_vaccinated_uids]) / rr(self.pars.dur_reach_peak / dt)
        self.ti_start_waning[new_vaccinated_uids] = sim.ti + rr(self.pars.dur_reach_peak / dt)

        # Update linear boost for HIV positives
        self.linear_boost[hiv_pos_uids] = self.linear_boost[hiv_pos_uids] * hiv_rel_linear_boost[np.digitize(hiv.cd4[hiv_pos_uids], cd4_bins)]
        duration_to_reach_peak = (self.peak_immunity[hiv_pos_uids] - self.immunity[hiv_pos_uids]) / self.linear_boost[hiv_pos_uids]
        self.ti_start_waning[hiv_pos_uids] = sim.ti + duration_to_reach_peak
        return

    def update_immunity_by_vaccination(self, sim):
        """
        Update Immunity levels for vaccinated agents
        1) Update immunity and transmission for all vaccinated agents: Distinguish by agents in their boosting and waning phases
        2) Update maternal transmission

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

            # Prevent infection and transmission params
            prevent_infection_param = self.pars.prevent_infection
            prevent_transmission_param = self.pars[f'prevent_transmission_{state}']

            if len(uids):
                ################################################################################
                # 1) Update immunity and transmission for all vaccinated agents

                # Uids to boost and uids to wane
                uids_to_wane = uids[self.ti_start_waning[uids] <= sim.ti]
                uids_to_boost = uids.remove(uids_to_wane)

                # Boosting
                if len(uids_to_boost):
                    # Set linear boost value, for HIV positives this will depend on CD4 counts
                    self.set_linear_boost(sim, uids_to_boost)

                    # Boost linearly to peak immunity
                    self.immunity[uids_to_boost] = np.minimum(self.peak_immunity[uids_to_boost], self.immunity[uids_to_boost] + self.linear_boost[uids_to_boost])
                    self.rel_sus_immunity[uids_to_boost] = prevent_infection_param * self.immunity[uids_to_boost]
                    self.rel_trans_immunity[uids_to_boost] = prevent_transmission_param * self.immunity[uids_to_boost]

                # Waning
                if len(uids_to_wane):
                    self.immunity[uids_to_wane] = self.logistic_decay(sim.ti - self.ti_start_waning[uids_to_wane], self.peak_immunity[uids_to_wane])
                    self.rel_sus_immunity[uids_to_wane] = prevent_infection_param * self.immunity[uids_to_wane]
                    self.rel_trans_immunity[uids_to_wane] = prevent_transmission_param * self.immunity[uids_to_wane]

                # Ensure values are non-negative
                self.rel_sus_immunity[uids] = np.minimum(1, self.rel_sus_immunity[uids]).clip(0)  # Make sure immunity is between 0 and 1
                self.rel_trans_immunity[uids] = np.minimum(1, self.rel_trans_immunity[uids]).clip(0)  # Make sure immunity is between 0 and 1

                ################################################################################
                # 2) Update transmission for maternal network
                self.rel_trans_immunity_maternal[uids] = self.rel_trans_immunity[uids] * (1 - self.pars.rel_efficacy_red_maternal)

                # Ensure values are non-negative
                self.rel_trans_immunity_maternal[uids] = np.minimum(1, self.rel_trans_immunity_maternal[uids]).clip(0)  # Make sure immunity is between 0 and 1

        return

    def compute_trans_sus(self, sim):
        syph = sim.diseases.syphilis
        rel_trans = syph.rel_trans * (1 - self.rel_trans_immunity)
        rel_trans_maternal = syph.rel_trans_maternal * (1 - self.rel_trans_immunity_maternal)
        rel_sus = syph.rel_sus * (1 - self.rel_sus_immunity)
        return rel_trans, rel_trans_maternal, rel_sus

    def update_dur_infection(self, sim):
        """
        Reduce ti_secondary, and ti_tertiary for vaccinated, infected agents.
        """
        # Extract parameters and indices
        syph = sim.diseases.syphilis
        has_primary_syphilis = syph.primary
        has_secondary_syphilis = syph.secondary
        ti_infected = syph.ti_infected
        is_vaccinated = self.vaccinated
        dur_inf_updated = self.dur_inf_updated
        # New reinfections
        new_reinfected = syph.reinfected & (syph.ti_infected == sim.ti - 1)

        # Get people that are syphilis infected and vaccinated, who hadn't had their infection duration updated
        # This ensures that vaccinated people who get infected get their duration updated as well as
        # infected people who get vaccinated.
        # Primary -> Secondary, Secondary -> Tertiary
        for state, next_state in zip(['primary', 'secondary'], ['secondary', 'latent']):
            ti_duration = getattr(syph, f'ti_{next_state}')
            # Get uids from agents whose duration haven't been updated ever, or who have been reinfected
            uids = (is_vaccinated & getattr(syph, state) & ti_duration.notnan & (~dur_inf_updated | new_reinfected)).uids

            current_duration = ti_duration[uids] - sim.ti
            reduce_duration_parameter = self.pars[f'reduce_dur_{state}']
            new_ti_duration = sim.ti + rr(current_duration * reduce_duration_parameter)
            ti_duration[uids] = new_ti_duration

            # Update values
            setattr(syph, f'ti_{next_state}', ti_duration)

            # Update bool - ensures we only update the infection durations once per agent
            self.dur_inf_updated[uids] = True

        return

    def update_latent_prognoses(self, sim):
        """
        A proportion of agents in the latent syphilis stage are scheduled to reactivate based on p_reactivation in the syphilis disease module.
        Here, we reduce the proportion of agents that are scheduled to reactivate to the secondary stage for vaccinated agents based
        on the reduce_p_reactivate parameter.
        Agents that have reactivation to secondary removed will stay latent for the rest of the simulation (unless treated)
        """
        syph = sim.diseases.syphilis
        vaccinated = self.vaccinated
        new_latent = (syph.ti_latent == sim.ti)

        # Scheduled to reactivate
        reactivate_bool = new_latent & (syph.ti_secondary >= sim.ti)
        reactive_uids = reactivate_bool.uids

        # Remove reactivation for a proportion of agents scheduled for reactivation according to reduce_p_reactivate param
        n_reactivate_vaccinated = len((reactivate_bool & vaccinated).uids)
        n_reactivate_vaccinated_target = syph.pars.p_reactivate.pars.p * self.pars.reduce_p_reactivate * len((vaccinated & new_latent).uids)
        n_remove_reactivation = n_reactivate_vaccinated - rr(n_reactivate_vaccinated_target)
        # Pick agents randomly
        if n_remove_reactivation > 0:
            bools = ss.random(strict=False).rvs(len((reactivate_bool & vaccinated).uids))
            choices = np.argsort(bools)[:n_remove_reactivation]
            target_uids = (reactivate_bool & vaccinated).uids[choices]

            # Reset ti_secondary
            syph.ti_secondary[target_uids] = np.nan

        return

    def update_peak_immunity(self, sim, uids):
        """
        Update peak immunity levels based on received doses
        """
        # Uids
        uids_first_dose = uids[self.doses[uids] == 1]
        uids_second_dose = uids[self.doses[uids] == 2]
        uids_third_dose = uids[self.doses[uids] == 3]

        # Peak immunities
        self.peak_immunity[uids_first_dose] = self.pars.first_dose_efficacy
        self.peak_immunity[uids_second_dose] = self.pars.second_dose_efficacy
        self.peak_immunity[uids_third_dose] = self.pars.third_dose_efficacy
        return

    def update_rel_sus_rel_trans(self, sim):
        """
        Update relative susceptibility and transmission, including transmission of maternal network
        """
        syph = sim.diseases.syphilis
        # Set rel trans and rel sus
        rel_trans, rel_trans_maternal, rel_sus = self.compute_trans_sus(sim)
        syph.rel_trans.set(uids=sim.people.auids, new_vals=rel_trans)
        syph.rel_trans_maternal.set(uids=sim.people.auids, new_vals=rel_trans_maternal)
        syph.rel_sus.set(uids=sim.people.auids, new_vals=rel_sus)

    def vaccinate(self, sim, uids, update_immunity_by_vaccination=True):
        """
        Vaccinate
        """
        # Set states
        self.vaccinated[uids] = True
        self.ti_vaccinated[uids] = sim.ti
        self.ti_nab_event[uids] = sim.ti

        # Reset start waning
        self.ti_start_waning[uids] = np.nan

        # Update number of doses
        self.doses[uids] += 1

        # Update peak immunity
        self.update_peak_immunity(sim, uids)

        # Update immunity
        if update_immunity_by_vaccination:
            # Update Immunity
            self.update_immunity_by_vaccination(sim)

        # Schedule a second dose for a proportion of agents who received their first vaccination
        dt = sim.dt_year
        first_dose_uids = uids[self.doses[uids] == 1]
        will_get_second_dose = first_dose_uids[self.pars.p_second_dose.rvs(first_dose_uids)]
        self.ti_second_dose[will_get_second_dose] = sim.ti + rr(self.pars.dose_interval.rvs(will_get_second_dose) / dt)

        # Schedule a third dose for a proportion of agents that have received their second dose
        second_dose_uids = uids[self.doses[uids] == 2]
        will_get_third_dose = second_dose_uids[self.pars.p_third_dose.rvs(second_dose_uids)]
        self.ti_third_dose[will_get_third_dose] = sim.ti + rr(self.pars.dose_interval.rvs(will_get_third_dose) / dt)

        return

    def update_results(self):
        """
        Update results
        """
        ti = self.sim.ti
        self.results['n_doses'][ti] = sum(self.doses.values)
        self.results['n_vaccinated'][ti] = len(self.vaccinated.uids)
        self.results['n_vaccinated_twice'][ti] = len((self.doses == 2).uids)
        self.results['n_vaccinated_triple'][ti] = len((self.doses == 3).uids)
        self.results['new_vaccinated'][ti] = len(((self.ti_vaccinated == ti) & (self.doses == 1)).uids)
        return

    def step(self):
        sim = self.sim
        syph = sim.diseases.syphilis
        # self.update_natural_immunity(sim) # Taking this out for now

        if np.floor(sim.timevec[sim.ti]) > self.start_year:
            # Get this time steps target coverage
            target_coverage = self.target_coverage[sim.ti]

            target_uids = self.get_targets(sim, self.num_doses, target_coverage)
            # If there are targets, vaccinate them and update immunity for all vaccinated agents
            if len(target_uids):
                self.vaccinate(sim, target_uids)
            # If there are no targets, only update immunity for all vaccinated agents
            else:
                self.update_immunity_by_vaccination(sim)

            # Update rel_sus, rel_trans and rel_trans_maternal
            self.update_rel_sus_rel_trans(sim)

            # Reduce duration of infection for vaccinated, infected agents
            self.update_dur_infection(sim)

            # Update latent prognoses by reducing p_reactivate for vaccinated, infected agents in latent state
            self.update_latent_prognoses(sim)

            # Update Results
            self.update_results()

        return
