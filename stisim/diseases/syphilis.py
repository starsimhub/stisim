"""
Define syphilis disease module
"""

import numpy as np
import sciris as sc
from sciris import randround as rr # Since used frequently
import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI

__all__ = ['Syphilis', 'SyphilisPlaceholder']


class SyphilisPlaceholder(ss.Disease):
    # A simple placeholder module to use when testing connectors

    def __init__(self, pars=None, **kwargs):
        super().__init__(name='syphilis')

        self.default_pars(
            prevalence=0.1,  # Target prevalance. If None, no automatic infections will be applied
        )
        self.update_pars(pars, **kwargs)
        self.add_states(
            ss.BoolArr('active'), # Active syphilis
            ss.FloatArr('ti_active'), # Time of active syphilis
        )
        self._prev_dist = ss.bernoulli(p=0)

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if not isinstance(self.pars.prevalence, sti.TimeSeries):
            ts = sti.TimeSeries(assumption=self.pars.prevalence)
        else:
            ts = self.pars.prevalence
        self._target_prevalence = ts.interpolate(sim.yearvec)

    def set_prognoses(self, target_uids, source_uids=None):
        self.active[target_uids] = True

    def update_pre(self):
        """
        When using a connector to the syphilis module, this is not needed. The connector should update the syphilis-positive state.
        """

        if self.pars.prevalence is None:
            return

        sim = self.sim

        # Get current prevalence
        n_active = self.active.count()
        prev = n_active/len(sim.people)
        target = self._target_prevalence[sim.ti]
        change = target-prev

        if change > 0:
            # Add a proportion of people that are not infected
            uids = self.active.false()
            self._prev_dist.set(p=change/(len(uids)/len(sim.people)))
            self.active[self._prev_dist.filter(uids)] = True
        elif change < 0:
            uids = self.active.true()
            self._prev_dist.set(p=-change/(len(uids)/len(sim.people)))
            self.active[self._prev_dist.filter(uids)] = False


class Syphilis(BaseSTI):

    def __init__(self, pars=None, init_prev_data=None, init_prev_latent_data=None, **kwargs):
        super().__init__()
        self.requires = 'structuredsexual'

        self.default_pars(
            # Adult syphilis natural history, all specified in years
            dur_primary = ss.uniform(low=3/52, high=10/52),  # https://pubmed.ncbi.nlm.nih.gov/9101629/
            dur_secondary = ss.lognorm_ex(mean=3.6/12, stdev=1.5/12),  # https://pubmed.ncbi.nlm.nih.gov/9101629/
            dur_early = ss.normal(1.5, 2/12),  # Assumption
            p_reactivate = ss.bernoulli(p=0.35),  # Probability of reactivating from latent to secondary
            time_to_reactivate = ss.lognorm_ex(mean=1, stdev=1),  # Time to reactivation
            p_tertiary = ss.bernoulli(p=0.35),  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917057/
            time_to_tertiary = ss.lognorm_ex(mean=20, stdev=8),  # Time to tertiary
            p_death = ss.bernoulli(p=0.05),  # probability of dying of tertiary syphilis
            time_to_death = ss.lognorm_ex(mean=5, stdev=5),  # Time to death

            # Transmission by stage
            beta=1.0,  # Placeholder
            beta_m2f=None,
            beta_f2m=None,
            beta_m2c=None,
            eff_condom=0.0,
            rel_trans_primary=1,
            rel_trans_secondary=1,
            rel_trans_latent=1,  # Baseline level; this decays exponentially with duration of latent infection
            rel_trans_tertiary=0.0,
            rel_trans_latent_half_life=1,

            # Congenital syphilis outcomes
            # Birth outcomes coded as:
            #   0: Miscarriage
            #   1: Neonatal death
            #   2: Stillborn
            #   3: Congenital syphilis
            #   4: Live birth without syphilis-related complications
            # Sources:
            #   - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5973824/)
            #   - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819963/
            birth_outcomes=sc.objdict(
                early = ss.choice(a=5, p=np.array([0.10, 0.15, 0.10, 0.50, 0.15])), # Outcomes for babies born to mothers with primary, secondary, or early latent infection
                late = ss.choice(a=5, p=np.array([0.05, 0.05, 0.05, 0.10, 0.75])), # Outcomes for babies born to mothers with late latent infection
            ),
            birth_outcome_keys=['miscarriage', 'nnd', 'stillborn', 'congenital'],
            anc_detection=0.8,

            # Initial conditions
            init_prev=ss.bernoulli(p=0),
            init_latent_prev=ss.bernoulli(p=0),
            dist_ti_init_infected=ss.uniform(low=-10 * 12, high=1),
            rel_init_prev=1,
        )

        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data
        self.init_prev_latent_data = init_prev_latent_data
        if init_prev_data is not None:
            self.pars.init_prev = ss.bernoulli(self.make_init_prev_fn)
        if init_prev_latent_data is not None:
            self.pars.init_latent_prev = ss.bernoulli(self.make_init_prev_latent_fn)

        self.add_states(
            # Adult syphilis states
            ss.BoolArr('primary'),      # Primary chancres
            ss.BoolArr('secondary'),    # Inclusive of those who may still have primary chancres
            ss.BoolArr('early'),        # Early latent
            ss.BoolArr('late'),         # Late latent
            ss.BoolArr('latent'),       # Can relapse to secondary, remain in latent, or progress to tertiary,
            ss.BoolArr('tertiary'),     # Includes complications (cardio/neuro/disfigurement)
            ss.BoolArr('immune'),       # After effective treatment people may acquire temp immunity
            ss.BoolArr('ever_exposed'), # Anyone ever exposed - stays true after treatment

            # Congenital syphilis states
            ss.BoolArr('congenital'),

            # Timestep of state changes
            ss.FloatArr('ti_transmitted'),
            ss.FloatArr('new_transmissions'),
            ss.FloatArr('cum_transmissions'),
            ss.FloatArr('ti_primary'),
            ss.FloatArr('ti_secondary'),
            ss.FloatArr('ti_latent'),
            ss.FloatArr('dur_early'),
            ss.FloatArr('ti_tertiary'),
            ss.FloatArr('ti_dead'),
            ss.FloatArr('ti_immune'),
            ss.FloatArr('ti_miscarriage'),
            ss.FloatArr('ti_nnd'),
            ss.FloatArr('ti_stillborn'),
            ss.FloatArr('ti_congenital'),
        )

        return

    @staticmethod
    def make_init_prev_fn(module, sim, uids):
        return sti.make_init_prev_fn(module, sim, uids, active=True)

    @staticmethod
    def make_init_prev_latent_fn(module, sim, uids):
        return sti.make_init_prev_fn(module, sim, uids, active=True, data=module.init_prev_latent_data)

    @property
    def naive(self):
        """ Never exposed """
        return ~self.ever_exposed

    @property
    def sus_not_naive(self):
        """ Susceptible but with syphilis antibodies, which persist after treatment """
        return self.susceptible & self.ever_exposed

    @property
    def active(self):
        """ Active infection includes primary and secondary stages """
        return self.primary | self.secondary

    @property
    def infectious(self):
        """ Infectious """
        return self.active | self.latent

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.pars.beta_m2f is not None:
            self.pars.beta['structuredsexual'][0] *= self.pars.beta_m2f
        if self.pars.beta_f2m is not None:
            self.pars.beta['structuredsexual'][1] *= self.pars.beta_f2m
        if self.pars.beta_m2c is not None:
            self.pars.beta['maternal'][0] *= self.pars.beta_m2c
        return

    def init_post(self):
        """ Make initial cases - TODO, figure out how to incorporate active syphilis here """
        initial_active_cases = self.pars.init_prev.filter()
        self.set_prognoses(initial_active_cases)
        still_sus = self.susceptible.uids

        # Natural history for initial latent cases
        initial_latent_cases = self.pars.init_latent_prev.filter(still_sus)
        ti_init_cases = self.pars.dist_ti_init_infected.rvs(initial_latent_cases).astype(int)
        self.set_prognoses(initial_latent_cases, ti=ti_init_cases)
        self.set_secondary_prognoses(initial_latent_cases)
        time_to_tertiary = self.pars.time_to_tertiary.rvs(initial_latent_cases)
        self.ti_tertiary[initial_latent_cases] = self.ti_latent[initial_latent_cases] + rr(time_to_tertiary / self.sim.dt)

        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        npts = self.sim.npts
        self.results += ss.Result(self.name, 'n_active', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'pregnant_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'detected_pregnant_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'adult_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'active_adult_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'active_prevalence', npts, dtype=float, scale=False)
        self.results += ss.Result(self.name, 'new_nnds', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_stillborns',  npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_congenital',  npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_congenital_deaths', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'cum_congenital',  npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'cum_congenital_deaths', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_deaths', npts, dtype=int, scale=True)

        # Add FSW and clients to results:
        self.results += ss.Result(self.name, 'prevalence_sw', npts, dtype=float)
        self.results += ss.Result(self.name, 'new_infections_sw', npts, dtype=float, scale=True)
        self.results += ss.Result(self.name, 'new_infections_not_sw', npts, dtype=float, scale=True)
        self.results += ss.Result(self.name, 'prevalence_client', npts, dtype=float)
        self.results += ss.Result(self.name, 'new_infections_client', npts, dtype=float, scale=True)
        self.results += ss.Result(self.name, 'new_infections_not_client', npts, dtype=float, scale=True)

        # Add risk groups to results
        for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
            for sex in ['female', 'male']:
                self.results += ss.Result(self.name, 'prevalence_risk_group_' + str(risk_group) + '_' + sex, npts, dtype=float)
                self.results += ss.Result(self.name, 'new_infections_risk_group_' + str(risk_group) + '_' + sex, npts, dtype=float, scale=True)

        # Add overall testing and treatment results, which might be assembled from numerous interventions
        self.results += ss.Result(self.name, 'new_false_pos', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_true_pos', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_false_neg', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_true_neg', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_treated_success', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_treated_failure', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_treated_unnecessary', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_treated', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_fetus_treated_success', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_fetus_treated_unnecessary', npts, dtype=int, scale=True)
        self.results += ss.Result(self.name, 'new_fetus_treated_failure', npts, dtype=int, scale=True)

        return

    def update_pre(self):
        """ Updates prior to interventions """
        ti = self.sim.ti
        dt = self.sim.dt

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Secondary from primary
        secondary_from_primary = self.primary & (self.ti_secondary <= ti)
        if len(secondary_from_primary.uids) > 0:
            self.secondary[secondary_from_primary] = True
            self.primary[secondary_from_primary] = False
            self.set_secondary_prognoses(secondary_from_primary.uids)

        # Hack to reset MultiRNGs in set_secondary_prognoses so they can be called again this timestep. TODO: Refactor
        self.pars.dur_secondary.jump(ti+1)

        # Secondary reactivation from latent
        secondary_from_latent = self.latent & (self.ti_latent >= ti) & (self.ti_secondary <= ti)
        if len(secondary_from_latent.uids) > 0:
            self.secondary[secondary_from_latent] = True
            self.latent[secondary_from_latent] = False
            self.set_secondary_prognoses(secondary_from_latent.uids)

        # Latent
        latent = self.secondary & (self.ti_latent <= ti)
        if len(latent.uids) > 0:
            self.latent[latent] = True
            self.secondary[latent] = False
            self.set_latent_prognoses(latent.uids)

        # Tertiary
        tertiary = self.latent & (self.ti_tertiary <= ti)
        self.tertiary[tertiary] = True
        self.latent[tertiary] = False

        # Trigger deaths
        deaths = (self.ti_dead <= ti).uids
        if len(deaths):
            self.sim.people.request_death(deaths)

        # Congenital syphilis deaths
        nnd = (self.ti_nnd <= ti).uids
        stillborn = (self.ti_stillborn <= ti).uids
        self.sim.people.request_death(nnd)
        self.sim.people.request_death(stillborn)

        # Congenital syphilis transmission outcomes
        congenital = (self.ti_congenital <= ti).uids
        self.congenital[congenital] = True
        self.susceptible[congenital] = False

        # Set rel_trans
        self.rel_trans[self.primary] = self.pars.rel_trans_primary
        self.rel_trans[self.secondary] = self.pars.rel_trans_secondary
        self.rel_trans[self.tertiary] = self.pars.rel_trans_tertiary
        # Latent rel_trans decays with duration of latent infection
        if len(self.latent.uids) > 0:
            self.set_latent_trans()

            # Set people to early
            lu = self.latent.uids
            is_early = lu[(ti - self.ti_latent[lu])*dt <= self.dur_early[lu]]
            is_late = lu[(ti - self.ti_latent[lu])*dt > self.dur_early[lu]]
            self.early[is_early] = True
            self.early[is_late] = False
            self.late[is_early] = False
            self.late[is_late] = True

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti

        self.results['n_active'][ti] = self.results['n_primary'][ti] + self.results['n_secondary'][ti]
        self.results['active_prevalence'][ti] = self.results['n_active'][ti] / np.count_nonzero(self.sim.people.alive)
        active_adults_num = len(((self.sim.people.age >= 15) & (self.sim.people.age < 50) & (self.active)).uids)
        infected_adults_num = len(((self.sim.people.age >= 15) & (self.sim.people.age < 50)& (self.infected)).uids)
        adults_denom = len(((self.sim.people.age >= 15) & (self.sim.people.age < 50)).uids)
        self.results['adult_prevalence'][ti] = infected_adults_num / adults_denom
        self.results['active_adult_prevalence'][ti] = active_adults_num / adults_denom

        # Pregnant women prevalence, if present
        if 'pregnancy' in self.sim.demographics.keys():
            preg_denom = np.count_nonzero(self.sim.people.pregnancy.pregnant)
            preg_num = np.count_nonzero(self.sim.people.pregnancy.pregnant & self.infected)
            detected_preg_num = preg_num*self.pars.anc_detection
            self.results['pregnant_prevalence'][ti] = sc.safedivide(preg_num, preg_denom)
            self.results['detected_pregnant_prevalence'][ti] = sc.safedivide(detected_preg_num, preg_denom)

        # Congenital results
        self.results['new_nnds'][ti]       = np.count_nonzero(self.ti_nnd == ti)
        self.results['new_stillborns'][ti] = np.count_nonzero(self.ti_stillborn == ti)
        self.results['new_congenital'][ti] = np.count_nonzero(self.ti_congenital == ti)
        self.results['new_congenital_deaths'][ti] = self.results['new_nnds'][ti] + self.results['new_stillborns'][ti]
        self.results['new_deaths'][ti] = np.count_nonzero(self.ti_dead == ti)

        # Add FSW and clients to results:
        fsw_infected = self.infected[self.sim.networks.structuredsexual.fsw]
        client_infected = self.infected[self.sim.networks.structuredsexual.client]
        if len(fsw_infected) > 0:
            self.results['prevalence_sw'][ti] = sum(fsw_infected) / len(fsw_infected)
            self.results['new_infections_sw'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.fsw).uids)
            self.results['new_infections_not_sw'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.fsw).uids)
        if len(client_infected) > 0:
            self.results['prevalence_client'][ti] = sum(client_infected) / len(client_infected)
            self.results['new_infections_client'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.client).uids)
            self.results['new_infections_not_client'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.client).uids)

        # Add risk groups
        for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
            for sex in ['female', 'male']:
                risk_group_infected = self.infected[(self.sim.networks.structuredsexual.risk_group == risk_group) & (self.sim.people[sex])]
                risk_group_new_inf = ((self.ti_infected == ti) & (self.sim.networks.structuredsexual.risk_group == risk_group) &  (self.sim.people[sex])).uids
                if len(risk_group_infected) > 0:
                    self.results['prevalence_risk_group_' + str(risk_group) + '_' + sex][ti] = sum(risk_group_infected) / len(risk_group_infected)
                    self.results['new_infections_risk_group_' + str(risk_group) + '_' + sex][ti] = len(risk_group_new_inf)

        return

    def finalize_results(self):
        self.results['cum_congenital'] = np.cumsum(self.results['new_congenital'])
        self.results['cum_congenital_deaths'] = np.cumsum(self.results['new_congenital_deaths'])
        super().finalize_results()
        return

    def set_latent_trans(self, ti=None):
        if ti is None: ti = self.sim.ti
        dt = self.sim.dt
        dur_latent = ti - self.ti_latent[self.latent]
        hl = self.pars.rel_trans_latent_half_life
        decay_rate = np.log(2) / hl if ~np.isnan(hl) else 0.
        latent_trans = self.pars.rel_trans_latent * np.exp(-decay_rate * dur_latent * dt)
        self.rel_trans[self.latent] = latent_trans
        return

    def make_new_cases(self):
        """
        Add new syphilis cases
        """
        targets, sources, networks = super().make_new_cases()
        for source, target, network in zip(sources, targets, networks):
            self.log.append(source, target, t=self.sim.year, network=self.sim.networks[network].name)
        return sources, targets, networks

    def set_prognoses(self, uids, source_uids=None, ti=None):
        """
        Set initial prognoses for adults newly infected with syphilis
        """
        if ti is None:
            ti = self.sim.ti
        else:
            # Check that ti is consistent with uids
            if not (sc.isnumber(ti) or len(ti) == len(uids)):
                errormsg = 'ti for set_prognoses must be int or array of length uids'
                raise ValueError(errormsg)

        dt = self.sim.dt
        self.new_transmissions[:] = 0  # Reset this every timestep

        # If someone has been infected by >1 person, remove duplicates
        uidx = np.unique(uids, return_index=True)[1]
        uids = ss.uids([uids[index] for index in sorted(uidx)])
        if source_uids is not None:
            source_uids = ss.uids([source_uids[index] for index in sorted(uidx)])
            unique_sources, counts = np.unique(source_uids, return_counts=True)
            self.ti_transmitted[unique_sources] = ti
            self.new_transmissions[unique_sources] = counts
            self.cum_transmissions[unique_sources] += counts

        self.susceptible[uids] = False
        self.ever_exposed[uids] = True
        self.primary[uids] = True
        self.infected[uids] = True
        self.ti_primary[uids] = ti
        self.ti_infected[uids] = ti

        # Primary to secondary
        dur_primary = self.pars.dur_primary.rvs(uids)
        self.ti_secondary[uids] = self.ti_primary[uids] + rr(dur_primary / dt)

        self.dur_early[uids] = self.pars.dur_early.rvs(uids)

        return

    def set_secondary_prognoses(self, uids):
        """ Set prognoses for people who have just progressed to secondary infection """
        dt = self.sim.dt
        dur_secondary = self.pars.dur_secondary.rvs(uids)
        self.ti_latent[uids] = self.ti_secondary[uids] + rr(dur_secondary / dt)
        return

    def set_latent_prognoses(self, uids):
        ti = self.sim.ti
        dt = self.sim.dt
        # Reactivators
        will_reactivate = self.pars.p_reactivate.rvs(uids)
        reactivate_uids = uids[will_reactivate]
        if len(reactivate_uids) > 0:
            time_to_reactivate = self.pars.time_to_reactivate.rvs(reactivate_uids)
            self.ti_secondary[reactivate_uids] = self.ti_latent[reactivate_uids] + rr(time_to_reactivate / dt)

        # Latent to tertiary
        nonreactivate_uids = uids[~will_reactivate]
        if len(nonreactivate_uids) > 0:
            is_tertiary = self.pars.p_tertiary.rvs(nonreactivate_uids)
            tertiary_uids = nonreactivate_uids[is_tertiary]
            if len(tertiary_uids) > 0:
                time_to_tertiary = self.pars.time_to_tertiary.rvs(tertiary_uids)
                self.ti_tertiary[tertiary_uids] = self.ti_latent[tertiary_uids] + rr(time_to_tertiary / dt)

                # Tertiary to dead
                will_die = self.pars.p_death.rvs(tertiary_uids)
                dead_uids = tertiary_uids[will_die]
                if len(dead_uids) > 0:
                    time_to_death = self.pars.time_to_death.rvs(dead_uids)
                    self.ti_dead[dead_uids] = self.ti_tertiary[dead_uids] + rr(time_to_death / dt)

        return

    def set_congenital(self, target_uids, source_uids=None):
        """
        Natural history of syphilis for congenital infection
        """
        ti = self.sim.ti
        dt = self.sim.dt

        # Determine outcomes
        for state in ['early', 'late']:

            source_state_inds = getattr(self, state)[source_uids].nonzero()[-1]
            uids = target_uids[source_state_inds]

            if len(uids) > 0:

                # Birth outcomes must be modified to add probability of susceptible birth
                birth_outcomes = self.pars.birth_outcomes[state]
                assigned_outcomes = birth_outcomes.rvs(uids)
                ages = self.sim.people.age

                # Schedule events
                for oi, outcome in enumerate(self.pars.birth_outcome_keys):
                    o_uids = uids[assigned_outcomes == oi]
                    if len(o_uids) > 0:
                        ti_outcome = f'ti_{outcome}'
                        vals = getattr(self, ti_outcome)
                        vals[o_uids] = ti + rr(-ages[o_uids] / dt)

                        setattr(self, ti_outcome, vals)

        return


