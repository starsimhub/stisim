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
            dur_exp2inf=ss.constant(1/52),

            # Symptoms and symptomatic testing
            p_symp=[
                ss.bernoulli(p=0.375),  # Women
                ss.bernoulli(p=0.375),  # Men
            ],
            p_symp_test=[
                ss.bernoulli(p=0.2),  # Women
                ss.bernoulli(p=0.01),  # Men
            ],
            dur_symp2test=[
                ss.lognorm_ex(1/12, 1/12),  # Women
                ss.lognorm_ex(1.5/12, 1/12),  # Men
            ],

            # PID and PID care-seeking
            p_pid=ss.bernoulli(p=0.2),
            dur_inf2pid=ss.lognorm_ex(1.5/12, 1/12),
            p_pid_test=ss.bernoulli(p=0.1),  # Women
            dur_pid2test=ss.lognorm_ex(0.5/12, 1/12),  # Women

            # Clearance
            dur_inf2clear=[
                ss.lognorm_ex(52/52, 5/52),  # Women
                ss.lognorm_ex(52/52, 5/52),  # Men
            ],

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
            ss.BoolArr('asymptomatic'),
            ss.BoolArr('symptomatic'),
            ss.BoolArr('pid'),
            ss.BoolArr('seeking_care'),
            ss.FloatArr('dur_inf'),
            ss.FloatArr('ti_exposed'),
            ss.FloatArr('ti_symptomatic'),
            ss.FloatArr('ti_seeks_care'),
            ss.FloatArr('ti_pid'),
            ss.FloatArr('ti_clearance'),
        )

        # Results by age and sex
        self.sex_results = sc.objdict()
        self.age_sex_results = sc.objdict()
        self.age_bins = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 100])
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
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        return

    def init_post(self):
        """ Make initial cases """
        super().init_post()
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        npts = self.sim.npts
        self.results += ss.Result(self.name, 'symp_prevalence', npts, dtype=float, scale=False, label="Symptomatic prevalence")
        self.results += ss.Result(self.name, 'incidence', npts, dtype=float, scale=False, label="Incidence")
        self.results += ss.Result(self.name, 'adult_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'symp_adult_prevalence', npts, dtype=float, scale=False)

        # Age/sex results
        for rkey in self.age_sex_result_keys:
            self.sex_results[rkey] = sc.objdict()
            self.age_sex_results[rkey] = sc.objdict()
            for skey in ['female', 'male', 'both']:
                self.sex_results[rkey][skey] = np.zeros(len(self.sim.yearvec))
                self.age_sex_results[rkey][skey] = np.zeros((len(self.age_bins)-1, len(self.sim.yearvec)))

        return
 
    def clear_infection(self, uids):
        self.exposed[uids] = False
        self.infected[uids] = False
        self.symptomatic[uids] = False
        self.asymptomatic[uids] = False
        self.pid[uids] = False
        self.seeking_care[uids] = False
        self.susceptible[uids] = True
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

        # Clear infections
        new_cleared = (self.infected & (self.ti_clearance <= ti)).uids
        self.clear_infection(new_cleared)

        # Progress symptoms
        new_symptomatic = (self.asymptomatic & (self.ti_symptomatic <= ti)).uids
        self.asymptomatic[new_symptomatic] = False
        self.symptomatic[new_symptomatic] = True
        self.ti_symptomatic[new_symptomatic] = ti

        # Progress PID
        new_pid = (~self.pid & (self.ti_pid <= ti)).uids
        self.pid[new_pid] = True
        self.ti_pid[new_pid] = ti

        # Symptomatic/PID care seeking
        old_seekers = (self.seeking_care).uids
        self.seeking_care[old_seekers] = False  # Remove the old
        self.ti_seeks_care[old_seekers] = np.nan  # Remove the old
        new_seekers = (~self.seeking_care & (self.ti_seeks_care <= ti)).uids
        self.seeking_care[new_seekers] = True
        self.ti_seeks_care[new_seekers] = ti

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        bins = self.age_bins

        self.results['symp_prevalence'][ti] = self.results['n_symptomatic'][ti] / np.count_nonzero(ppl.alive)
        self.results['incidence'][ti] = self.results['new_infections'][ti] / self.results['n_susceptible'][ti]
        adults = (self.sim.people.age >= 15) & (self.sim.people.age < 50)
        infected_adults = adults & self.infected
        symptomatic_adults = adults & self.symptomatic
        self.results['adult_prevalence'][ti] = np.count_nonzero(infected_adults) / np.count_nonzero(adults)
        self.results['symp_adult_prevalence'][ti] = np.count_nonzero(symptomatic_adults) / np.count_nonzero(adults)

        rmap = {'alive': 'both', 'female': 'female', 'male': 'male'}

        # Incidence and prevalence by age and sex
        for pkey, rkey in rmap.items():
            new_inf = ((self.ti_infected == ti) & ppl[pkey]).uids
            new_inf_ages = ppl.age[new_inf]
            n_sus = (self.susceptible & ppl[pkey]).uids
            n_sus_ages = ppl.age[n_sus]
            num, _ = np.histogram(new_inf_ages, bins=bins)
            denom, _ = np.histogram(n_sus_ages, bins=bins)
            self.age_sex_results['incidence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['incidence'][rkey][ti] = sc.safedivide(len(new_inf), len(n_sus))

        # Prevalence by age and sex
        for pkey, rkey in rmap.items():
            n_inf = (self.infected & ppl[pkey]).uids
            n_inf_ages = ppl.age[n_inf]
            num, _ = np.histogram(n_inf_ages, bins=bins)
            denom, _ = np.histogram(ppl.age[ppl[pkey]], bins=bins)
            self.age_sex_results['prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['prevalence'][rkey][ti] = sc.safedivide(len(n_inf), np.count_nonzero(ppl[pkey]))

            n_symp = (self.symptomatic & ppl[pkey]).uids
            n_symp_ages = ppl.age[n_symp]
            num, _ = np.histogram(n_symp_ages, bins=bins)
            self.age_sex_results['symp_prevalence'][rkey][:, ti] = sc.safedivide(num, denom)
            self.sex_results['symp_prevalence'][rkey][ti] = sc.safedivide(len(n_symp), np.count_nonzero(ppl[pkey]))

        return

    def finalize_results(self):
        super().finalize_results()
        return

    def set_exposure(self, uids):
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.asymptomatic[uids] = True
        self.ti_exposed[uids] = self.sim.ti
        dur_exp = self.pars.dur_exp2inf.rvs(uids)
        self.ti_infected[uids] = self.sim.ti + dur_exp/self.sim.dt
        return

    def set_symptoms(self, p, f_uids, m_uids):
        f_symp = p.p_symp[0].filter(f_uids)
        m_symp = p.p_symp[1].filter(m_uids)
        self.ti_symptomatic[f_symp] = self.ti_infected[f_symp]
        self.ti_symptomatic[m_symp] = self.ti_infected[m_symp]
        return f_symp, m_symp

    def set_duration(self, p, f_uids, m_uids):
        dur_inf_f = p.dur_inf2clear[0].rvs(f_uids)
        dur_inf_m = p.dur_inf2clear[1].rvs(m_uids)
        self.dur_inf[f_uids] = dur_inf_f
        self.dur_inf[m_uids] = dur_inf_m
        return

    def set_care_seeking(self, p, f_symp, m_symp):
        f_symp_test = p.p_symp_test[0].filter(f_symp)
        m_symp_test = p.p_symp_test[1].filter(m_symp)
        f_dur_symp2test = np.minimum(p.dur_symp2test[0].rvs(f_symp_test), self.dur_inf[f_symp_test])
        m_dur_symp2test = np.minimum(p.dur_symp2test[1].rvs(m_symp_test), self.dur_inf[m_symp_test])
        self.ti_seeks_care[f_symp_test] = self.ti_infected[f_symp_test] + f_dur_symp2test/self.sim.dt
        self.ti_seeks_care[m_symp_test] = self.ti_infected[m_symp_test] + m_dur_symp2test/self.sim.dt
        return

    def set_pid(self, p, f_uids):
        pid = p.p_pid.filter(f_uids)
        dur_prepid = np.minimum(p.dur_inf2pid.rvs(pid), self.dur_inf[pid])
        self.ti_pid[pid] = self.ti_infected[pid] + dur_prepid/self.sim.dt
        return pid

    def set_pid_care_seeking(self, p, pid):
        dt = self.sim.dt
        pid_test = p.p_pid_test.filter(pid)
        dur_pid2test = np.minimum(p.dur_pid2test.rvs(pid_test), self.dur_inf[pid_test])
        self.ti_seeks_care[pid_test] = np.minimum(self.ti_seeks_care[pid_test], self.ti_infected[pid_test] + dur_pid2test/dt)
        return

    def wipe_dates(self, uids):
        """ Clear all previous dates """
        self.ti_exposed[uids] = np.nan
        self.ti_infected[uids] = np.nan
        self.ti_symptomatic[uids] = np.nan
        self.ti_pid[uids] = np.nan
        self.ti_seeks_care[uids] = np.nan
        self.ti_clearance[uids] = np.nan
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

        self.set_exposure(uids)
        f_symp, m_symp = self.set_symptoms(p, f_uids, m_uids)
        self.set_duration(p, f_uids, m_uids)
        self.set_care_seeking(p, f_symp, m_symp)
        pid = self.set_pid(p, f_uids)
        self.set_pid_care_seeking(p, pid)

        # Determine when people recover
        self.ti_clearance[uids] = self.ti_infected[uids] + self.dur_inf[uids]/self.sim.dt

        return

