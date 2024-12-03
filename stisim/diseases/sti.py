"""
Template class for SEIS-type STIs
Used for chlamydia, gonorrhea, and trich
"""

import numpy as np
import starsim as ss
import sciris as sc
from stisim.utils import make_init_prev_fn
ss_int_ = ss.dtypes.int

__all__ = ['BaseSTI', 'SEIS']


def true(arr):
    return np.count_nonzero(arr)


class BaseSTI(ss.Infection):
    """
    Base class for sexually transmitted infections.
    Modifies make_new_cases to account for barrier protection.
    """
    def __init__(self, name=None, pars=None, init_prev_data=None, **kwargs):
        super().__init__(name=name)
        self.requires = 'structuredsexual'
        self.define_pars(
            unit='month',
            beta=0,  # Placeholder: no transmission. This will be set in validate_beta
            eff_condom=1,
            rel_init_prev=1,
        )
        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data
        if init_prev_data is not None:
            self.pars.init_prev = ss.bernoulli(self.make_init_prev_fn)

        return

    @staticmethod
    def make_init_prev_fn(module, sim, uids):
        return make_init_prev_fn(module, sim, uids, active=True)

    def validate_beta(self, run_checks=False):
        betamap = super().validate_beta(run_checks=run_checks)
        if self.pars.beta_m2f is not None and betamap and 'structuredsexual' in betamap.keys():
            betamap['structuredsexual'][0] = self.pars.beta_m2f
            betamap['structuredsexual'][1] = self.pars.beta_m2f * self.pars.rel_beta_f2m
        if self.pars.beta_m2c is not None and betamap and 'maternal' in betamap.keys():
            betamap['maternal'][0] = ss.beta(self.pars.beta_m2c, 'month').init(parent=self.sim.t)
        return betamap

    def infect(self):
        """ Determine who gets infected on this timestep via transmission on the network """

        new_cases = []
        sources = []
        networks = []
        betamap = self.validate_beta()

        rel_trans = self.rel_trans.asnew(self.infectious * self.rel_trans)
        rel_sus   = self.rel_sus.asnew(self.susceptible * self.rel_sus)

        for i, (nkey,net) in enumerate(self.sim.networks.items()):
            nk = ss.standardize_netkey(nkey)
            if len(net): # Skip networks with no edges
                edges = net.edges
                p1p2b0 = [edges.p1, edges.p2, betamap[nk][0]] # Person 1, person 2, beta 0
                p2p1b1 = [edges.p2, edges.p1, betamap[nk][1]] # Person 2, person 1, beta 1
                for src, trg, beta in [p1p2b0, p2p1b1]:
                    if beta: # Skip networks with no transmission
                        beta_per_dt = net.net_beta(disease_beta=beta, disease=self) # Compute beta for this network and timestep
                        randvals = self.trans_rng.rvs(src, trg) # Generate a new random number based on the two other random numbers
                        args = (src, trg, rel_trans, rel_sus, beta_per_dt, randvals) # Set up the arguments to calculate transmission
                        target_uids, source_uids = self.compute_transmission(*args) # Actually calculate it
                        new_cases.append(target_uids)
                        sources.append(source_uids)
                        networks.append(np.full(len(target_uids), dtype=ss_int_, fill_value=i))

        # Finalize
        if len(new_cases) and len(sources):
            new_cases = ss.uids.cat(new_cases)
            new_cases, inds = new_cases.unique(return_index=True)
            sources = ss.uids.cat(sources)[inds]
            networks = np.concatenate(networks)[inds]
        else:
            new_cases = ss.uids()
            sources = ss.uids()
            networks = np.empty(0, dtype=ss_int_)

        return new_cases, sources, networks



class SEIS(BaseSTI):

    def __init__(self, pars=None, name=None, init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)
        self.define_pars(
            # Settings
            include_care=True,  # Determines whether testing results are included

            # Natural history
            dur_exp=ss.constant(ss.dur(1, 'week')),  # Initial latent period: how long after exposure before you can infect others

            # Symptoms and symptomatic testing
            p_symp=[
                ss.bernoulli(p=0.375),  # Women
                ss.bernoulli(p=0.375),  # Men
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(12, 'week')),   # Women
                ss.lognorm_ex(ss.dur(0.25, 'week'), ss.dur(3, 'week')), # Men
            ],
            p_symp_clear=[
                ss.bernoulli(p=0.0),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            dur_symp=[
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(26, 'week')),  # Duration of symptoms
                ss.lognorm_ex(ss.dur(1, 'week'), ss.dur(26, 'week')),  # Duration of symptoms
            ],
            p_symp_care=[
                ss.bernoulli(p=0.3),  # Reset for each individual STI
                ss.bernoulli(p=0.2),  # As above
            ],
            dur_symp2care=[  # For those who test, how long before they seek care - reset for each individual STI
                ss.lognorm_ex(ss.dur(4, 'week'), ss.dur(4, 'week')),  # Women
                ss.lognorm_ex(ss.dur(6, 'week'), ss.dur(4, 'week')),  # Men
            ],

            # PID and PID care-seeking
            p_pid=ss.bernoulli(p=0.2),
            dur_prepid=ss.lognorm_ex(ss.dur(6, 'week'), ss.dur(4, 'week')),
            p_pid_care=ss.bernoulli(p=0.1),  # Women
            dur_pid2care=ss.lognorm_ex(ss.dur(2, 'week'), ss.dur(4, 'week')),  # Women

            # Clearance
            dur_asymp2clear=[  # Duration of untreated asymptomatic infection (excl initial latent)
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Women
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Men
            ],
            dur_symp2clear=[  # Duration of untreated symptomatic infection (excl initial latent)
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Women
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Men
            ],
            dur_postsymp2clear=[
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Women
                ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),  # Men
            ],
            dur_pid2clear=ss.lognorm_ex(ss.dur(52, 'week'), ss.dur(5, 'week')),

            # Transmission. In the parent class, beta is set to 1. Here, we set beta_m2f and beta_m2c
            beta_m2f=None,
            rel_beta_f2m=0.5,
            beta_m2c=None,

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        self.define_states(
            # Natural history
            ss.State('exposed'),
            ss.State('infected'),
            ss.State('asymptomatic'),
            ss.State('symptomatic'),
            ss.State('postsymptomatic'),
            ss.State('pid'),
            ss.State('seeking_care'),
            ss.FloatArr('dur_inf'),
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_symptomatic'),
            ss.FloatArr('ti_symp_clear'),
            ss.FloatArr('ti_seeks_care'),
            ss.FloatArr('ti_pid'),
            ss.FloatArr('ti_clearance'),
        )

        # Results by age and sex
        self.age_sex_results = sc.objdict()
        self.age_bins = np.array([0, 15, 20, 25, 35, 65, 100])
        self.age_sex_result_keys = [
            'incidence',
            'prevalence',
            'symp_prevalence',
        ]

        return

    @property
    def treatable(self):
        """ Active bacterial presence -- includes exposed and infected, and responds to treatment """
        return self.exposed | self.infected

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = [
            ss.Result('female_prevalence', scale=False, label="Prevalence - F"),
            ss.Result('male_prevalence', scale=False, label="Prevalence - M"),

            ss.Result('symp_prevalence', scale=False, label="Symptomatic prevalence"),
            ss.Result('female_symp_prevalence', scale=False, label="Symptomatic prevalence - F"),
            ss.Result('male_symp_prevalence', scale=False, label="Symptomatic prevalence - M"),

            ss.Result('adult_prevalence', scale=False, label="Adult prevalence"),
            ss.Result('female_adult_prevalence', scale=False, label="Adult prevalence - F"),
            ss.Result('male_adult_prevalence', scale=False, label="Adult prevalence - F"),

            ss.Result('incidence', scale=False, label="Incidence"),
            ss.Result('female_incidence', scale=False, label="Incidence - F"),
            ss.Result('male_incidence', scale=False, label="Incidence - M"),

            ss.Result('new_female_infections', dtype=int, label="New infections - F"),
            ss.Result('new_male_infections', dtype=int, label="New infections - F"),

            ss.Result('symp_adult_prevalence', scale=False, label="Symptomatic adult prevalence"),
            ss.Result('female_symp_adult_prevalence', scale=False, label="Symptomatic adult female prevalence"),
            ss.Result('male_symp_adult_prevalence', scale=False, label="Symptomatic adult female prevalence"),

            ss.Result('n_female_infected', dtype=int, label="Number infected - F"),
            ss.Result('n_male_infected', dtype=int, label="Number infected - F"),
            ss.Result('n_female_symptomatic', dtype=int, label="Number symptomatic - F"),
            ss.Result('n_male_symptomatic', dtype=int, label="Number symptomatic - F"),

            ss.Result('new_symptomatic', dtype=int, label="New symptomatic"),
            ss.Result('new_female_symptomatic', dtype=int, label="New symptomatic - F"),
            ss.Result('new_male_symptomatic', dtype=int, label="New symptomatic - F"),
        ]

        if self.pars.include_care:
            results += [
                ss.Result('new_care_seekers', dtype=int, label="New care seekers"),

                # Add overall testing and treatment results, which might be assembled from numerous interventions
                ss.Result('new_false_pos', dtype=int, label="New false positives"),
                ss.Result('new_true_pos', dtype=int, label="New true positives"),
                ss.Result('new_false_neg', dtype=int, label="New false negatives"),
                ss.Result('new_true_neg', dtype=int, label="New true negatives"),
                ss.Result('new_false_pos_f', dtype=int, label="New false positives - F"),
                ss.Result('new_true_pos_f', dtype=int, label="New true positives - F"),
                ss.Result('new_false_neg_f', dtype=int, label="New false negatives - F"),
                ss.Result('new_true_neg_f', dtype=int, label="New true negatives - F"),
                ss.Result('new_false_pos_m', dtype=int, label="New false positives - M"),
                ss.Result('new_true_pos_m', dtype=int, label="New true positives - M"),
                ss.Result('new_false_neg_m', dtype=int, label="New false negatives - M"),
                ss.Result('new_true_neg_m', dtype=int, label="New true negatives - M"),
                ss.Result('new_treated_success', dtype=int, label="Successful treatments"),
                ss.Result('new_treated_failure', dtype=int, label="Unsuccessful treatments"),
                ss.Result('new_treated_unnecessary', dtype=int, label="Unnecessary treatments"),
                ss.Result('new_treated_success_symp', dtype=int, label="Successful treatments (symptomatic)"),
                ss.Result('new_treated_success_asymp', dtype=int, label="Successful treatments (asymptomatic)"),
                ss.Result('new_treated', dtype=int, label="Treatments"),
            ]

        self.define_results(*results)

        # Age/sex results
        for rkey in self.age_sex_result_keys:
            # self.sex_results[rkey] = sc.objdict()
            self.age_sex_results[rkey] = sc.objdict()
            for skey in ['female', 'male', 'both']:
                # self.sex_results[rkey][skey] = np.zeros(len(self.sim.timevec))
                self.age_sex_results[rkey][skey] = np.zeros((len(self.age_bins)-1, self.t.npts))

        return

    def clear_infection(self, uids):
        self.exposed[uids] = False
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.postsymptomatic[uids] = False
        self.pid[uids] = False
        self.seeking_care[uids] = False
        self.susceptible[uids] = True
        past_care_seekers = uids[(self.ti_seeks_care[uids] < self.sim.ti).nonzero()[-1]]
        self.ti_seeks_care[past_care_seekers] = np.nan
        self.ti_clearance[uids] = self.sim.ti

    def step_state(self):
        """ Updates for this timestep """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Exposed -> infected
        new_infected = (self.exposed & (self.ti_infected <= ti)).uids
        if len(new_infected):
            self.exposed[new_infected] = False
            self.infected[new_infected] = True
            self.ti_infected[new_infected] = ti

        # Presymptomatic -> symptomatic
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        if len(new_symptomatic):
            self.asymptomatic[new_symptomatic] = False
            self.symptomatic[new_symptomatic] = True
            self.ti_symptomatic[new_symptomatic] = ti

        # Symptomatic -> post-symptomatic
        new_postsymptomatic = (self.symptomatic & (self.ti_symp_clear <= ti)).uids
        if len(new_postsymptomatic):
            self.symptomatic[new_postsymptomatic] = False
            self.postsymptomatic[new_postsymptomatic] = True
            self.ti_symp_clear[new_postsymptomatic] = ti

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Progress PID
        new_pid = (~self.pid & (self.ti_pid <= ti)).uids
        self.pid[new_pid] = True
        self.ti_pid[new_pid] = ti

        # Symptomatic/PID care seeking
        old_seekers = (self.seeking_care).uids
        self.seeking_care[old_seekers] = False
        self.ti_seeks_care[old_seekers] = np.nan  # Remove the old
        new_seekers = (self.infected & (self.ti_seeks_care <= ti)).uids
        self.seeking_care[new_seekers] = True
        self.ti_seeks_care[new_seekers] = ti
        # If we don't update this results here, it's possible that someone could seek care,
        # then get reinfected, and the reinfection would wipe all their dates so we'd miss
        # counting them in the results.
        if self.pars.include_care:
            self.results['new_care_seekers'][ti] = true(self.ti_seeks_care == ti)

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        bins = self.age_bins

        adults = (self.sim.people.age >= 15) & (self.sim.people.age <= 65)
        women = adults & self.sim.people.female
        men = adults & self.sim.people.male
        infected_adults = adults & self.infected
        infected_women = women & self.infected
        infected_men = men & self.infected

        self.results['female_prevalence'][ti] = true(self.infected & self.sim.people.female) / true(self.sim.people.female)
        self.results['male_prevalence'][ti] = true(self.infected & self.sim.people.male) / true(self.sim.people.male)

        self.results['symp_prevalence'][ti] = true(self.symptomatic) / true(self.sim.people.alive)
        self.results['female_symp_prevalence'][ti] = true(self.symptomatic & self.sim.people.female) / true(self.sim.people.female)
        self.results['male_symp_prevalence'][ti] = true(self.symptomatic & self.sim.people.male) / true(self.sim.people.male)

        self.results['adult_prevalence'][ti] = true(infected_adults) / true(adults)
        self.results['female_adult_prevalence'][ti] = true(infected_women) / true(women)
        self.results['male_adult_prevalence'][ti] = true(infected_men) / true(men)

        self.results['new_female_infections'][ti] = true((self.ti_infected == ti) & self.sim.people.female)
        self.results['new_male_infections'][ti] = true((self.ti_infected == ti) & self.sim.people.male)

        self.results['symp_adult_prevalence'][ti] = true(self.symptomatic & adults) / true(adults)
        self.results['female_symp_adult_prevalence'][ti] = true(self.symptomatic & women) / true(women)
        self.results['male_symp_adult_prevalence'][ti] = true(self.symptomatic & men) / true(men)

        self.results['n_female_infected'][ti] = true(self.infected & women)
        self.results['n_male_infected'][ti] = true(self.infected & men)
        self.results['n_female_symptomatic'][ti] = true(self.symptomatic & women)
        self.results['n_male_symptomatic'][ti] = true(self.symptomatic & men)

        self.results['new_symptomatic'][ti] = true(self.ti_symptomatic == ti)
        self.results['new_female_symptomatic'][ti] = true((self.ti_symptomatic == ti) & self.sim.people.female)
        self.results['new_male_symptomatic'][ti] = true((self.ti_symptomatic == ti) & self.sim.people.male)

        self.results['incidence'][ti] = sc.safedivide(self.results['new_infections'][ti], self.results['n_susceptible'][ti])

        # rmap = {'alive': 'both', 'female': 'female', 'male': 'male'}
        rmap = {'female': 'female', 'male': 'male'}

        # Incidence and prevalence by age and sex
        for pkey, rkey in rmap.items():
            new_inf = ((self.ti_infected == ti) & ppl[pkey]).uids
            new_inf_ages = ppl.age[new_inf]
            n_sus = (self.susceptible & ppl[pkey]).uids
            n_sus_ages = ppl.age[n_sus]
            num, _ = np.histogram(new_inf_ages, bins=bins)
            denom, _ = np.histogram(n_sus_ages, bins=bins)
            self.age_sex_results['incidence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.results[rkey+'_incidence'][ti] = sc.safedivide(len(new_inf), len(n_sus))

        # Prevalence by age and sex
        for pkey, rkey in rmap.items():
            n_inf = (self.infected & ppl[pkey]).uids
            n_inf_ages = ppl.age[n_inf]
            num, _ = np.histogram(n_inf_ages, bins=bins)
            denom, _ = np.histogram(ppl.age[ppl[pkey]], bins=bins)
            self.age_sex_results['prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.results[rkey+'_prevalence'][ti] = sc.safedivide(len(n_inf), np.count_nonzero(ppl[pkey]))

            n_symp = (self.symptomatic & ppl[pkey]).uids
            n_symp_ages = ppl.age[n_symp]
            num, _ = np.histogram(n_symp_ages, bins=bins)
            self.age_sex_results['symp_prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.results[rkey+'_symp_prevalence'][ti] = sc.safedivide(len(n_symp), np.count_nonzero(ppl[pkey]))

        return

    def set_exposure(self, uids):
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.asymptomatic[uids] = True
        self.ti_exposed[uids] = self.sim.ti
        dur_exp = self.pars.dur_exp.rvs(uids)
        self.ti_infected[uids] = self.sim.ti + dur_exp
        return

    def set_symptoms(self, p, f_uids, m_uids):
        f_symp, f_asymp = p.p_symp[0].split(f_uids)
        m_symp, m_asymp = p.p_symp[1].split(m_uids)
        f_dur_presymp = self.pars.dur_presymp[0].rvs(f_symp)
        m_dur_presymp = self.pars.dur_presymp[1].rvs(m_symp)
        self.ti_symptomatic[f_symp] = self.ti_infected[f_symp] + f_dur_presymp
        self.ti_symptomatic[m_symp] = self.ti_infected[m_symp] + m_dur_presymp
        return f_symp, m_symp, f_asymp, m_asymp

    def set_symp_clearance(self, p, f_symp, m_symp):
        f_symp_clear, f_symp_persist = p.p_symp_clear[0].split(f_symp)
        m_symp_clear, m_symp_persist = p.p_symp_clear[1].split(m_symp)
        f_dur_symp = self.pars.dur_symp[0].rvs(f_symp_clear)
        m_dur_symp = self.pars.dur_symp[1].rvs(m_symp_clear)
        self.ti_symp_clear[f_symp_clear] = self.ti_symptomatic[f_symp_clear] + f_dur_symp
        self.ti_symp_clear[m_symp_clear] = self.ti_symptomatic[m_symp_clear] + m_dur_symp
        return f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist

    def set_care_seeking(self, p, f_symp, m_symp):
        f_symp_care = p.p_symp_care[0].filter(f_symp)
        m_symp_care = p.p_symp_care[1].filter(m_symp)
        f_dur_symp2care = p.dur_symp2care[0].rvs(f_symp_care)
        m_dur_symp2care = p.dur_symp2care[1].rvs(m_symp_care)
        self.ti_seeks_care[f_symp_care] = self.ti_symptomatic[f_symp_care] + f_dur_symp2care
        self.ti_seeks_care[m_symp_care] = self.ti_symptomatic[m_symp_care] + m_dur_symp2care
        return

    def set_pid(self, p, f_uids):
        pid = p.p_pid.filter(f_uids)
        dur_prepid = p.dur_prepid.rvs(pid)
        self.ti_pid[pid] = self.ti_infected[pid] + dur_prepid
        return pid

    def set_pid_care_seeking(self, p, pid):
        pid_care = p.p_pid_care.filter(pid)
        dur_pid2care = p.dur_pid2care.rvs(pid_care)
        self.ti_seeks_care[pid_care] = np.minimum(self.ti_seeks_care[pid_care], self.ti_infected[pid_care] + dur_pid2care)
        return

    def set_duration(self, p, f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist, f_asymp, m_asymp, pid):

        # Duration of infection for those with persistant symptoms, transient symptoms, and asymptomatic infection
        dur_inf_f_symp_clear = p.dur_postsymp2clear[0].rvs(f_symp_clear)
        dur_inf_m_symp_clear = p.dur_postsymp2clear[1].rvs(m_symp_clear)
        dur_inf_f_symp_persist = p.dur_symp2clear[0].rvs(f_symp_persist)
        dur_inf_m_symp_persist = p.dur_symp2clear[1].rvs(m_symp_persist)
        dur_inf_f_asymp = p.dur_asymp2clear[0].rvs(f_asymp)
        dur_inf_m_asymp = p.dur_asymp2clear[1].rvs(m_asymp)
        dur_inf_pid = p.dur_pid2clear.rvs(pid)
        self.ti_clearance[f_symp_clear] = dur_inf_f_symp_clear + self.ti_symp_clear[f_symp_clear]
        self.ti_clearance[m_symp_clear] = dur_inf_m_symp_clear + self.ti_symp_clear[m_symp_clear]
        self.ti_clearance[f_symp_persist] = dur_inf_f_symp_persist + self.ti_symptomatic[f_symp_persist]
        self.ti_clearance[m_symp_persist] = dur_inf_m_symp_persist + self.ti_symptomatic[m_symp_persist]
        self.ti_clearance[f_asymp] = dur_inf_f_asymp + self.ti_infected[f_asymp]
        self.ti_clearance[m_asymp] = dur_inf_m_asymp + self.ti_infected[m_asymp]
        self.ti_clearance[pid] = np.maximum(self.ti_clearance[pid], dur_inf_pid + self.ti_pid[pid])
        return

    def wipe_dates(self, uids):
        """ Clear all previous dates """
        self.ti_exposed[uids] = np.nan
        self.ti_infected[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        self.ti_symp_clear[uids] = np.nan
        self.ti_pid[uids] = np.nan
        # self.ti_seeks_care[uids] = np.nan
        self.ti_clearance[uids] = np.nan
        self.dur_inf[uids] = np.nan
        return

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected
        """
        super().set_prognoses(uids, source_uids)
        self.wipe_dates(uids)

        ppl = self.sim.people
        p = self.pars
        f_uids = ppl.female.uids.intersect(uids)
        m_uids = ppl.male.uids.intersect(uids)

        self.set_exposure(uids)  # Set exposure
        f_symp, m_symp, f_asymp, m_asymp = self.set_symptoms(p, f_uids, m_uids)  # Set symptoms & presymptomatic duration
        f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist = self.set_symp_clearance(p, f_symp, m_symp)
        self.set_care_seeking(p, f_symp, m_symp)  # Determine who seeks care and when
        pid = self.set_pid(p, f_uids)  # Determine who developes PID and when
        self.set_pid_care_seeking(p, pid)
        self.set_duration(p, f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist, f_asymp, m_asymp, pid)

        # Determine overall duration of infection
        self.dur_inf[uids] = self.ti_clearance[uids] - self.ti_infected[uids]

        if (self.dur_inf[uids] < 0).any():
            errormsg = 'Invalid durations of infection'
            raise ValueError(errormsg)

        return

