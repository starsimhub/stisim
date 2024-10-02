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


class BaseSTI(ss.Infection):
    """
    Base class for sexually transmitted infections.
    Modifies make_new_cases to account for barrier protection.
    """
    def __init__(self, pars=None, init_prev_data=None, **kwargs):
        super().__init__()
        self.requires = 'structuredsexual'
        self.default_pars(
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

    def make_new_cases(self):
        """
        Create new cases via contact networks. Most of this is copied from the Starsim class,
        but the main difference is that beta_per_dt takes the disease module.
        """
        new_cases = []
        sources = []
        networks = []
        betamap = self._check_betas()

        for i, (nkey,net) in enumerate(self.sim.networks.items()):
            if not len(net):
                break

            nbetas = betamap[nkey]
            edges = net.edges

            rel_trans = self.rel_trans.asnew(self.infectious * self.rel_trans)
            rel_sus   = self.rel_sus.asnew(self.susceptible * self.rel_sus)
            p1p2b0 = [edges.p1, edges.p2, nbetas[0]]
            p2p1b1 = [edges.p2, edges.p1, nbetas[1]]
            for src, trg, beta in [p1p2b0, p2p1b1]:

                # Skip networks with no transmission
                if beta == 0:
                    continue

                # Calculate probability of a->b transmission.
                if net.postnatal or net.prenatal:
                    beta_per_dt = net.beta_per_dt(disease_beta=beta, dt=self.sim.dt)
                else:
                    beta_per_dt = net.beta_per_dt(disease_beta=beta, dt=self.sim.dt, disease=self)
                p_transmit = rel_trans[src] * rel_sus[trg] * beta_per_dt

                # Generate a new random number based on the two other random numbers
                rvs_s = self.rng_source.rvs(src)
                rvs_t = self.rng_target.rvs(trg)
                rvs = ss.combine_rands(rvs_s, rvs_t)

                new_cases_bool = rvs < p_transmit
                new_cases.append(trg[new_cases_bool])
                sources.append(src[new_cases_bool])
                networks.append(np.full(np.count_nonzero(new_cases_bool), dtype=ss_int_, fill_value=i))

        # Tidy up
        if len(new_cases) and len(sources):
            new_cases = ss.uids.cat(new_cases)
            new_cases, inds = new_cases.unique(return_index=True)
            sources = ss.uids.cat(sources)[inds]
            networks = np.concatenate(networks)[inds]
        else:
            new_cases = np.empty(0, dtype=int)
            sources = np.empty(0, dtype=int)
            networks = np.empty(0, dtype=int)

        if len(new_cases):
            self._set_cases(new_cases, sources)

        return new_cases, sources, networks


class SEIS(BaseSTI):

    def __init__(self, pars=None, name=None, init_prev_data=None, **kwargs):
        super().__init__(name=name, init_prev_data=init_prev_data)
        self.default_pars(

            # Natural history
            dur_exp=ss.constant(1/52),  # Initial latent period: how long after exposure before you can infect others

            # Symptoms and symptomatic testing
            p_symp=[
                ss.bernoulli(p=0.375),  # Women
                ss.bernoulli(p=0.375),  # Men
            ],
            dur_presymp=[  # For those who develop symptoms, how long before symptoms appear
                ss.lognorm_ex(1/52, 12/52),  # Women:
                ss.lognorm_ex(0.5/52, 12/52),  # Men: symptoms should appear within days
            ],
            p_symp_clear=[
                ss.bernoulli(p=0.0),  # Women
                ss.bernoulli(p=0.0),  # Men
            ],
            dur_symp=[
                ss.lognorm_ex(1/52, 26/52),  # Duration of symptoms
                ss.lognorm_ex(1/52, 26/52),  # Duration of symptoms
            ],
            p_symp_care=[
                ss.bernoulli(p=0.3),  # Women: 0.2*1+0.8*.2
                ss.bernoulli(p=0.2),  # Men: 0.5*1+0.5*.25
            ],
            dur_symp2care=[  # For those who test, how long before they seek care
                ss.lognorm_ex(1/12, 1/12),  # Women
                ss.lognorm_ex(1.5/12, 1/12),  # Men
            ],

            # PID and PID care-seeking
            p_pid=ss.bernoulli(p=0.2),
            dur_prepid=ss.lognorm_ex(1.5/12, 1/12),
            p_pid_care=ss.bernoulli(p=0.1),  # Women
            dur_pid2care=ss.lognorm_ex(0.5/12, 1/12),  # Women

            # Clearance
            dur_asymp2clear=[  # Duration of untreated asymptomatic infection (excl initial latent)
                ss.lognorm_ex(52/52, 5/52),  # Women
                ss.lognorm_ex(52/52, 5/52),  # Men
            ],
            dur_symp2clear=[  # Duration of untreated symptomatic infection (excl initial latent)
                ss.lognorm_ex(52/52, 5/52),  # Women
                ss.lognorm_ex(52/52, 5/52),  # Men
            ],
            dur_postsymp2clear=[
                ss.lognorm_ex(52/52, 5/52),  # Women
                ss.lognorm_ex(52/52, 5/52),  # Men
            ],
            dur_pid2clear=ss.lognorm_ex(52/52, 5/52),

            # Transmission
            beta=1.0,  # Placeholder
            beta_m2f=None,
            beta_f2m=None,
            beta_m2c=None,

            # Initial conditions
            init_prev=ss.bernoulli(p=0.01)
        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            # Natural history
            ss.BoolArr('exposed'),
            ss.BoolArr('infected'),
            ss.BoolArr('asymptomatic'),
            ss.BoolArr('symptomatic'),
            ss.BoolArr('postsymptomatic'),
            ss.BoolArr('pid'),
            ss.BoolArr('seeking_care'),
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

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.beta_m2f is not None:
            self.pars.beta['structuredsexual'][0] *= self.pars.beta_m2f
            if self.pars.beta_f2m is None:
                self.pars.beta_f2m = self.pars.beta_m2f / 2
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        if self.pars.beta_m2c is not None:
            self.pars.beta['maternalnet'][1] *= self.pars.beta_m2c
        return

    def init_post(self):
        """ Make initial cases """
        super().init_post()
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        npts = self.sim.npts
        self.results += ss.Result(self.name, 'female_prevalence', npts, dtype=float, scale=False, label="Prevalence - F")
        self.results += ss.Result(self.name, 'male_prevalence', npts, dtype=float, scale=False, label="Prevalence - M")

        self.results += ss.Result(self.name, 'symp_prevalence', npts, dtype=float, scale=False, label="Symptomatic prevalence")
        self.results += ss.Result(self.name, 'female_symp_prevalence', npts, dtype=float, scale=False, label="Symptomatic prevalence - F")
        self.results += ss.Result(self.name, 'male_symp_prevalence', npts, dtype=float, scale=False, label="Symptomatic prevalence - M")

        self.results += ss.Result(self.name, 'adult_prevalence', npts, dtype=float, scale=False, label="Adult prevalence")
        self.results += ss.Result(self.name, 'female_adult_prevalence', npts, dtype=float, scale=False, label="Adult prevalence - F")
        self.results += ss.Result(self.name, 'male_adult_prevalence', npts, dtype=float, scale=False, label="Adult prevalence - F")

        self.results += ss.Result(self.name, 'incidence', npts, dtype=float, scale=False, label="Incidence")
        self.results += ss.Result(self.name, 'female_incidence', npts, dtype=float, scale=False, label="Incidence - F")
        self.results += ss.Result(self.name, 'male_incidence', npts, dtype=float, scale=False, label="Incidence - M")

        self.results += ss.Result(self.name, 'new_female_infections', npts, dtype=float, scale=False, label="New infections - F")
        self.results += ss.Result(self.name, 'new_male_infections', npts, dtype=float, scale=False, label="New infections - F")

        self.results += ss.Result(self.name, 'symp_adult_prevalence', npts, dtype=float, scale=False, label="Symptomatic adult prevalence")
        self.results += ss.Result(self.name, 'female_symp_adult_prevalence', npts, dtype=float, scale=False, label="Symptomatic adult female prevalence")
        self.results += ss.Result(self.name, 'male_symp_adult_prevalence', npts, dtype=float, scale=False, label="Symptomatic adult female prevalence")

        self.results += ss.Result(self.name, 'n_female_infected', npts, dtype=int, scale=True, label="Number infected - F")
        self.results += ss.Result(self.name, 'n_male_infected', npts, dtype=int, scale=True, label="Number infected - F")
        self.results += ss.Result(self.name, 'n_female_symptomatic', npts, dtype=int, scale=True, label="Number symptomatic - F")
        self.results += ss.Result(self.name, 'n_male_symptomatic', npts, dtype=int, scale=True, label="Number symptomatic - F")

        self.results += ss.Result(self.name, 'new_symptomatic', npts, dtype=int, scale=True, label="New symptomatic")
        self.results += ss.Result(self.name, 'new_female_symptomatic', npts, dtype=int, scale=True, label="New symptomatic - F")
        self.results += ss.Result(self.name, 'new_male_symptomatic', npts, dtype=int, scale=True, label="New symptomatic - F")

        self.results += ss.Result(self.name, 'new_care_seekers', npts, dtype=int, scale=True, label="New care seekers")

        # Add overall testing and treatment results, which might be assembled from numerous interventions
        self.results += ss.Result(self.name, 'new_false_pos', npts, dtype=int, scale=True, label="New false positives")
        self.results += ss.Result(self.name, 'new_true_pos', npts, dtype=int, scale=True, label="New true positives")
        self.results += ss.Result(self.name, 'new_false_neg', npts, dtype=int, scale=True, label="New false negatives")
        self.results += ss.Result(self.name, 'new_true_neg', npts, dtype=int, scale=True, label="New true negatives")
        self.results += ss.Result(self.name, 'new_false_pos_f', npts, dtype=int, scale=True, label="New false positives - F")
        self.results += ss.Result(self.name, 'new_true_pos_f', npts, dtype=int, scale=True, label="New true positives - F")
        self.results += ss.Result(self.name, 'new_false_neg_f', npts, dtype=int, scale=True, label="New false negatives - F")
        self.results += ss.Result(self.name, 'new_true_neg_f', npts, dtype=int, scale=True, label="New true negatives - F")
        self.results += ss.Result(self.name, 'new_false_pos_m', npts, dtype=int, scale=True, label="New false positives - M")
        self.results += ss.Result(self.name, 'new_true_pos_m', npts, dtype=int, scale=True, label="New true positives - M")
        self.results += ss.Result(self.name, 'new_false_neg_m', npts, dtype=int, scale=True, label="New false negatives - M")
        self.results += ss.Result(self.name, 'new_true_neg_m', npts, dtype=int, scale=True, label="New true negatives - M")
        self.results += ss.Result(self.name, 'new_treated_success', npts, dtype=int, scale=True, label="Successful treatments")
        self.results += ss.Result(self.name, 'new_treated_failure', npts, dtype=int, scale=True, label="Unsuccessful treatments")
        self.results += ss.Result(self.name, 'new_treated_unnecessary', npts, dtype=int, scale=True, label="Unnecessary treatments")
        self.results += ss.Result(self.name, 'new_treated_success_symp', npts, dtype=int, scale=True, label="Successful treatments (symptomatic)")
        self.results += ss.Result(self.name, 'new_treated_success_asymp', npts, dtype=int, scale=True, label="Successful treatments (asymptomatic)")
        self.results += ss.Result(self.name, 'new_treated', npts, dtype=int, scale=True, label="Treatments")

        # Age/sex results
        for rkey in self.age_sex_result_keys:
            # self.sex_results[rkey] = sc.objdict()
            self.age_sex_results[rkey] = sc.objdict()
            for skey in ['female', 'male', 'both']:
                # self.sex_results[rkey][skey] = np.zeros(len(self.sim.yearvec))
                self.age_sex_results[rkey][skey] = np.zeros((len(self.age_bins)-1, len(self.sim.yearvec)))

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
        past_care_seekers = uids[(self.ti_seeks_care[uids]<self.sim.ti).nonzero()[-1]]
        self.ti_seeks_care[past_care_seekers] = np.nan
        self.ti_clearance[uids] = self.sim.ti

    def update_pre(self):
        """ Updates prior to interventions """
        ti = self.sim.ti

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
        self.results['new_care_seekers'][ti] = np.count_nonzero(self.ti_seeks_care == ti)

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

        self.results['female_prevalence'][ti] = np.count_nonzero(self.infected & self.sim.people.female) / np.count_nonzero(self.sim.people.female)
        self.results['male_prevalence'][ti] = np.count_nonzero(self.infected & self.sim.people.male) / np.count_nonzero(self.sim.people.male)

        self.results['symp_prevalence'][ti] = np.count_nonzero(self.symptomatic) / np.count_nonzero(self.sim.people.alive)
        self.results['female_symp_prevalence'][ti] = np.count_nonzero(self.symptomatic & self.sim.people.female) / np.count_nonzero(self.sim.people.female)
        self.results['male_symp_prevalence'][ti] = np.count_nonzero(self.symptomatic & self.sim.people.male) / np.count_nonzero(self.sim.people.male)

        self.results['adult_prevalence'][ti] = np.count_nonzero(self.infected & adults) / np.count_nonzero(adults)
        self.results['female_adult_prevalence'][ti] = np.count_nonzero(self.infected & women) / np.count_nonzero(women)
        self.results['male_adult_prevalence'][ti] = np.count_nonzero(self.infected & men) / np.count_nonzero(men)

        self.results['new_female_infections'][ti] = np.count_nonzero((self.ti_infected == ti) & self.sim.people.female)
        self.results['new_male_infections'][ti] = np.count_nonzero((self.ti_infected == ti) & self.sim.people.male)

        self.results['symp_adult_prevalence'][ti] = np.count_nonzero(self.symptomatic & adults) / np.count_nonzero(adults)
        self.results['female_symp_adult_prevalence'][ti] = np.count_nonzero(self.symptomatic & women) / np.count_nonzero(women)
        self.results['male_symp_adult_prevalence'][ti] = np.count_nonzero(self.symptomatic & men) / np.count_nonzero(men)

        self.results['n_female_infected'][ti] = np.count_nonzero(self.infected & women)
        self.results['n_male_infected'][ti] = np.count_nonzero(self.infected & men)
        self.results['n_female_symptomatic'][ti] = np.count_nonzero(self.symptomatic & women)
        self.results['n_male_symptomatic'][ti] = np.count_nonzero(self.symptomatic & men)

        self.results['new_symptomatic'][ti] = np.count_nonzero(self.ti_symptomatic == ti)
        self.results['new_female_symptomatic'][ti] = np.count_nonzero((self.ti_symptomatic == ti) & self.sim.people.female)
        self.results['new_male_symptomatic'][ti] = np.count_nonzero((self.ti_symptomatic == ti) & self.sim.people.male)

        self.results['incidence'][ti] = self.results['new_infections'][ti] / self.results['n_susceptible'][ti]

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

    def finalize_results(self):
        super().finalize_results()
        return

    def set_exposure(self, uids):
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.asymptomatic[uids] = True
        self.ti_exposed[uids] = self.sim.ti
        dur_exp = self.pars.dur_exp.rvs(uids)
        self.ti_infected[uids] = self.sim.ti + dur_exp/self.sim.dt
        return

    def set_symptoms(self, p, f_uids, m_uids):
        f_symp, f_asymp = p.p_symp[0].split(f_uids)
        m_symp, m_asymp = p.p_symp[1].split(m_uids)
        f_dur_presymp = self.pars.dur_presymp[0].rvs(f_symp)
        m_dur_presymp = self.pars.dur_presymp[1].rvs(m_symp)
        self.ti_symptomatic[f_symp] = self.ti_infected[f_symp] + f_dur_presymp/self.sim.dt
        self.ti_symptomatic[m_symp] = self.ti_infected[m_symp] + m_dur_presymp/self.sim.dt
        return f_symp, m_symp, f_asymp, m_asymp

    def set_symp_clearance(self, p, f_symp, m_symp):
        f_symp_clear, f_symp_persist = p.p_symp_clear[0].split(f_symp)
        m_symp_clear, m_symp_persist = p.p_symp_clear[1].split(m_symp)
        f_dur_symp = self.pars.dur_symp[0].rvs(f_symp_clear)
        m_dur_symp = self.pars.dur_symp[1].rvs(m_symp_clear)
        self.ti_symp_clear[f_symp_clear] = self.ti_symptomatic[f_symp_clear] + f_dur_symp/self.sim.dt
        self.ti_symp_clear[m_symp_clear] = self.ti_symptomatic[m_symp_clear] + m_dur_symp/self.sim.dt
        return f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist

    def set_care_seeking(self, p, f_symp, m_symp):
        f_symp_care = p.p_symp_care[0].filter(f_symp)
        m_symp_care = p.p_symp_care[1].filter(m_symp)
        f_dur_symp2care = p.dur_symp2care[0].rvs(f_symp_care)
        m_dur_symp2care = p.dur_symp2care[1].rvs(m_symp_care)
        self.ti_seeks_care[f_symp_care] = self.ti_symptomatic[f_symp_care] + f_dur_symp2care/self.sim.dt
        self.ti_seeks_care[m_symp_care] = self.ti_symptomatic[m_symp_care] + m_dur_symp2care/self.sim.dt
        return

    def set_pid(self, p, f_uids):
        pid = p.p_pid.filter(f_uids)
        dur_prepid = p.dur_prepid.rvs(pid)
        self.ti_pid[pid] = self.ti_infected[pid] + dur_prepid/self.sim.dt
        return pid

    def set_pid_care_seeking(self, p, pid):
        dt = self.sim.dt
        pid_care = p.p_pid_care.filter(pid)
        dur_pid2care = p.dur_pid2care.rvs(pid_care)
        self.ti_seeks_care[pid_care] = np.minimum(self.ti_seeks_care[pid_care], self.ti_infected[pid_care] + dur_pid2care/dt)
        return

    def set_duration(self, p, f_symp_clear, m_symp_clear, f_symp_persist, m_symp_persist, f_asymp, m_asymp, pid):
        dt = self.sim.dt

        # Duration of infection for those with persistant symptoms, transient symptoms, and asymptomatic infection
        dur_inf_f_symp_clear = p.dur_postsymp2clear[0].rvs(f_symp_clear)
        dur_inf_m_symp_clear = p.dur_postsymp2clear[1].rvs(m_symp_clear)
        dur_inf_f_symp_persist = p.dur_symp2clear[0].rvs(f_symp_persist)
        dur_inf_m_symp_persist = p.dur_symp2clear[1].rvs(m_symp_persist)
        dur_inf_f_asymp = p.dur_asymp2clear[0].rvs(f_asymp)
        dur_inf_m_asymp = p.dur_asymp2clear[1].rvs(m_asymp)
        dur_inf_pid = p.dur_pid2clear.rvs(pid)
        self.ti_clearance[f_symp_clear] = dur_inf_f_symp_clear/dt + self.ti_symp_clear[f_symp_clear]
        self.ti_clearance[m_symp_clear] = dur_inf_m_symp_clear/dt + self.ti_symp_clear[m_symp_clear]
        self.ti_clearance[f_symp_persist] = dur_inf_f_symp_persist/dt + self.ti_symptomatic[f_symp_persist]
        self.ti_clearance[m_symp_persist] = dur_inf_m_symp_persist/dt + self.ti_symptomatic[m_symp_persist]
        self.ti_clearance[f_asymp] = dur_inf_f_asymp/dt + self.ti_infected[f_asymp]
        self.ti_clearance[m_asymp] = dur_inf_m_asymp/dt + self.ti_infected[m_asymp]
        self.ti_clearance[pid] = np.maximum(self.ti_clearance[pid], dur_inf_pid/dt + self.ti_pid[pid])
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
        self.dur_inf[uids] = (self.ti_clearance[uids] - self.ti_infected[uids])*self.sim.dt

        if (self.dur_inf[uids] < 0).any():
            errormsg = 'Invalid durations of infection'
            raise ValueError(errormsg)

        return

