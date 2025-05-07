"""
Define syphilis interventions for STIsim
"""

import starsim as ss
import numpy as np
import pandas as pd
import sciris as sc
from stisim.interventions.base_interventions import STIDx, STITest, STITreatment
from stisim.utils import TimeSeries

__all__ = ["SyphDx", "SyphTx", "NewbornTreatment", "SyphTest", "ANCSyphTest", "NewbornSyphTest"]


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
            treat_eff=ss.bernoulli(p=.98),
            fetus_treat_eff=ss.bernoulli(p=.98),
            fetus_age_cutoff_treat_eff=-0.25,  # Reduced treatment efficacy for fetuses in the last trimester
            treat_eff_reduced=ss.bernoulli(p=.75)  # Reduced efficacy for fetuses older than cut off
        )
        self.update_pars(pars, **kwargs)
        return

    def change_states(self, disease, treat_succ):
        """ Change the states of people who are treated """
        self.sim.diseases[disease].primary[treat_succ] = False
        self.sim.diseases[disease].secondary[treat_succ] = False
        self.sim.diseases[disease].early[treat_succ] = False
        self.sim.diseases[disease].late[treat_succ] = False
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
        treat_st_cutoff = self.pars.fetus_treat_eff.filter(fetus_st_cutoff)

        # Combine
        treat_uids = treat_gt_cutoff.concat(treat_st_cutoff)

        # Store results - Failure
        sim.diseases.syphilis.results['new_fetus_treated_failure'][self.ti] = len(fetus_gt_cutoff) + len(fetus_st_cutoff) - len(treat_uids)
        # Store results - Unnecessary. Subtract successful uids below
        sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][self.ti] += len(treat_uids)

        # Change birth outcomes for successfully treated unborn babies
        # If fetus is not yet infected, this will do nothing
        if len(treat_uids):
            # Birth outcomes must be modified to add probability of susceptible birth
            birth_outcomes = sim.diseases.syphilis.pars.birth_outcome_keys

            # Change states
            for oi, outcome in enumerate(birth_outcomes):
                if outcome != 'normal':
                    ti_outcome = f'ti_{outcome}'
                    vals = getattr(sim.diseases.syphilis, ti_outcome)
                    successful_uids = treat_uids & vals.notnan.uids
                    vals[successful_uids] = np.nan
                    setattr(sim.diseases.syphilis, ti_outcome, vals)

                # Store results - Success
                sim.diseases.syphilis.results['new_fetus_treated_success'][self.ti] += len(successful_uids)
                # Store results - Unnecessary. Subtract the uids that were successful
                sim.diseases.syphilis.results['new_fetus_treated_unnecessary'][self.ti] -= len(successful_uids)

        return

    def step(self):
        """
        Apply treatment
        """
        sim = self.sim
        treat_uids = super().step()
        # Treat unborn babies of successfully treated mothers
        treat_pregnant_uids = sim.people.pregnancy.pregnant.uids & self.outcomes['successful']
        overtreat_pregnant_uids = sim.people.pregnancy.pregnant.uids & self.outcomes['unnecessary']
        # Store results - Unnecessary
        sim.diseases.syphilis.results['new_treated_unnecessary_pregnant'][self.ti] += len(overtreat_pregnant_uids)
        if len(treat_pregnant_uids):
            self.treat_fetus(sim, mother_uids=treat_pregnant_uids)
        return treat_uids


class NewbornTreatment(SyphTx):

    def __init__(self, pars=None, eligibility=None, max_capacity=None, years=None, name=None, *args, **kwargs):
        super().__init__(name=name, eligibility=eligibility, max_capacity=max_capacity, years=years, *args)
        self.define_pars(
            by_sex=False,
        )
        self.update_pars(pars, **kwargs)
        return

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
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None, name=None, label=None, newborn_test=None, **kwargs):
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, eligibility=eligibility, product=product, name=name, label=label)
        self.define_pars(
            dt_scale=True,
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
                dd[group] = sc.smoothinterp(self.t.yearvec, list(vals.keys()), list(vals.values()), smoothness=0)
            return dd
        else: return self.test_prob_data

    @staticmethod
    def make_test_prob_fn(self, sim, uids):
        """ Process symptomatic testing probabilites over time by sex and risk group """

        if sc.isnumber(self.test_prob_data):
            test_prob = self.test_prob_data
        elif isinstance(self.test_prob_data, TimeSeries):
            test_prob = self.test_prob_data.interpolate(sim.now)
        elif sc.checktype(self.test_prob_data, 'arraylike'):
            # year_ind = sc.findnearest(self.years, self.t.now('year'))
            # self.test_prob_data[year_ind]
            all_vals = sc.smoothinterp(self.t.yearvec, self.years, self.test_prob_data, smoothness=0)
            test_prob = all_vals[self.ti]
        elif isinstance(self.test_prob_data, dict):
            test_prob = pd.Series(index=uids)
            n_risk_groups = sim.networks.structuredsexual.pars.n_risk_groups
            for rg in range(n_risk_groups):
                for sex in ['female', 'male']:
                    for sw in [0, 1]:
                        dkey = (rg, sex, sw)
                        conditions = (sim.people[sex] & (sim.networks.structuredsexual.risk_group==rg))
                        if sw:
                            if sex == 'female': conditions = conditions & sim.networks.structuredsexual.fsw
                            if sex == 'male':   conditions = conditions & sim.networks.structuredsexual.client
                        test_prob[conditions[uids]] = self.test_prob_data[dkey][self.ti]
        else:
            errormsg = 'Format of test_prob_data must be float, array, or dict.'
            raise ValueError(errormsg)

        # Scale and validate
        test_prob = test_prob * self.pars.rel_test
        if self.pars.dt_scale:
            test_prob = test_prob * self.t.dt_year
        test_prob = np.clip(test_prob, a_min=0, a_max=1)

        return test_prob

    def step(self, uids=None):
        super().step(uids=uids)
        sim = self.sim

        after_start = (self.t.now('year') >= self.start)
        before_stop = (self.t.now('year') < self.stop)

        if (after_start & before_stop):
            # Schedule newborn tests if the mother is positive
            if self.newborn_test is not None:
                new_pos = self.ti_positive == self.ti
                if new_pos.any():
                    pos_mother_inds = np.in1d(sim.networks.maternalnet.p1, new_pos.uids)
                    unborn_uids = sim.networks.maternalnet.p2[pos_mother_inds]
                    ti_births = sim.networks.maternalnet.edges.end[pos_mother_inds].astype(int)
                    self.newborn_test.schedule(unborn_uids, ti_births)

        return

    def update_results(self):
        ti = self.ti
        super().update_results()
        if 'positive' in self.outcomes.keys() and (self.ti_positive == ti).any():
            self.update_positives()
        if 'negative' in self.outcomes.keys() and (self.ti_negative == ti).any():
            self.update_negatives()

        return

    def update_positives(self):
        ti = self.ti
        new_pos = self.ti_positive == ti
        # Count true/false positives
        false_pos = np.count_nonzero(self.sim.diseases.syphilis.susceptible[new_pos])
        true_pos = np.count_nonzero(self.sim.diseases.syphilis.infected[new_pos])
        self.sim.diseases.syphilis.results['new_false_pos'][ti] += false_pos
        self.sim.diseases.syphilis.results['new_true_pos'][ti] += true_pos
        return

    def update_negatives(self):
        ti = self.ti
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
    def __init__(self, test_prob_data=None, years=None, start=None, stop=None, pars=None, product=None, eligibility=None, name=None, label=None, newborn_test=None, **kwargs):
        super().__init__(test_prob_data=test_prob_data, years=years, start=start, stop=stop, eligibility=eligibility, product=product, name=name, label=label, **kwargs)
        self.define_pars(
            dt_scale=False,
        )
        self.update_pars(pars, **kwargs)
        self.test_timing = ss.randint(1, 8)
        if self.eligibility is None:
            self.eligibility = lambda sim: sim.demographics.pregnancy.pregnant
        return

    def get_testers(self, sim):
        # For ANC testing, only administer scheduled tests
        return (self.ti_scheduled == self.ti).uids

    def schedule_tests(self):
        """ Schedule a test for newly pregnant women """
        sim = self.sim
        newly_preg = (sim.demographics.pregnancy.ti_pregnant == (self.ti)).uids
        self.test_prob.pars['p'] = self.make_test_prob_fn(self, sim, newly_preg)
        will_test = self.test_prob.filter(newly_preg)
        ti_test = self.ti + self.test_timing.rvs(will_test)
        self.ti_scheduled[will_test] = ti_test

    def step(self):
        self.schedule_tests()  # Check for newly pregnant women so they can be added to the schedule
        super().step()
        return


class NewbornSyphTest(SyphTest):
    """
    Test given to newborns if the mother was confirmed to have syphilis at any stage of the pregnancy
    """

    def get_testers(self, sim):
        # For newborn testing, only administer scheduled tests, and account for probability of testing at this point
        eligible_uids = (self.ti_scheduled == self.ti).uids
        accept_uids = self.test_prob.filter(eligible_uids)
        return accept_uids

