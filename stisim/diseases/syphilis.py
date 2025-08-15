"""
Define syphilis disease module
"""

import numpy as np
import sciris as sc
from sciris import randround as rr # Since used frequently
import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI, BaseSTIPars

__all__ = ['Syphilis', 'SyphilisPlaceholder']


# Define some helper functions
def count(arr): return np.count_nonzero(arr)
def div(a, b): return sc.safedivide(a, b)
def countdiv(a, b): return sc.safedivide(count(a), count(b))
def cond_prob(a, b): return sc.safedivide(count(a & b), count(b))


class SyphilisPlaceholder(ss.Disease):
    # A simple placeholder module to use when testing connectors

    def __init__(self, pars=None, **kwargs):
        super().__init__(name='syphilis')

        self.define_pars(
            prevalence=0.1,  # Target prevalance. If None, no automatic infections will be applied
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
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
        self._target_prevalence = ts.interpolate(sim.timevec)

    def set_prognoses(self, target_uids, source_uids=None):
        self.active[target_uids] = True

    def step_state(self):
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


class SyphPars(BaseSTIPars):
    def __init__(self, **kwargs):
        # Adult syphilis natural history
        self.dur_primary = ss.normal(ss.weeks(6), ss.weeks(1))  # https://pubmed.ncbi.nlm.nih.gov/9101629/
        self.dur_secondary = ss.lognorm_ex(ss.months(3.6), ss.months(1.5))  # https://pubmed.ncbi.nlm.nih.gov/9101629/
        self.dur_early = ss.uniform(ss.months(12), ss.months(14))  # Assumption
        self.p_reactivate = ss.bernoulli(p=0.35)  # Probability of reactivating from latent to secondary
        self.time_to_reactivate = ss.lognorm_ex(ss.years(1), ss.years(1))  # Time to reactivation
        self.p_tertiary = ss.bernoulli(p=0.35)  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917057/
        self.time_to_tertiary = ss.normal(ss.years(20), ss.years(2))  # Time to tertiary
        self.p_death = ss.bernoulli(p=0.05)  # probability of dying of tertiary syphilis
        self.time_to_death = ss.lognorm_ex(ss.years(5), ss.years(5))  # Time to death

        # Transmission by stage
        self.eff_condom = 0.0
        self.rel_trans_primary = 1
        self.rel_trans_secondary = 1
        self.rel_trans_latent = 1  # Baseline level; this decays exponentially with duration of latent infection
        self.rel_trans_tertiary = 0.0
        self.rel_trans_latent_half_life = ss.years(1)

        # Congenital syphilis outcomes
        # Birth outcomes coded as:
        #   0: Miscarriage
        #   1: Neonatal death
        #   2: Stillborn
        #   3: Congenital syphilis
        #   4: Live birth without syphilis-related complications - may be preterm or low birth weight
        # Sources:
        #   - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5973824/)
        #   - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819963/
        self.birth_outcomes = sc.objdict(
            active=ss.choice(a=5, p=np.array([0.00, 0.10, 0.20, 0.45, 0.25])),  # Outcomes for babies born to mothers with primary or secondary infection
            early=ss.choice(a=5, p=np.array([0.00, 0.05, 0.10, 0.40, 0.45])),  # Outcomes for babies born to mothers with early latent infection
            late=ss.choice(a=5, p=np.array([0.00, 0.00, 0.10, 0.10, 0.80])),  # Outcomes for babies born to mothers with late latent infection
        )
        self.birth_outcome_keys = ['miscarriage', 'nnd', 'stillborn', 'congenital', 'normal']
        self.anc_detection = 0.8

        # Initial conditions
        self.init_prev = ss.bernoulli(p=0)
        self.init_latent_prev = ss.bernoulli(p=0)
        self.dist_ti_init_infected = ss.constant(0)  # Experimented with negative values, but safer to use 0
        self.rel_init_prev = 1

        # Update
        self.update(kwargs)


class Syphilis(BaseSTI):

    def __init__(self, pars=None, name='syph', init_prev_data=None, init_prev_latent_data=None, **kwargs):
        super().__init__(name=name)

        # Define default parameters
        default_pars = SyphPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data
        self.init_prev_latent_data = init_prev_latent_data
        if init_prev_data is not None:
            self.pars.init_prev = ss.bernoulli(self.make_init_prev_fn)
        if init_prev_latent_data is not None:
            self.pars.init_latent_prev = ss.bernoulli(self.make_init_prev_latent_fn)

        # Whether to store detailed results
        self.store_sw = False
        self.store_risk_groups = False

        self.define_states(
            # Adult syphilis states
            ss.BoolState('primary'),      # Primary chancres
            ss.BoolState('secondary'),    # Inclusive of those who may still have primary chancres
            ss.BoolState('early'),        # Early latent
            ss.BoolState('late'),         # Late latent
            ss.BoolState('latent'),       # Can relapse to secondary, remain in latent, or progress to tertiary,
            ss.BoolState('tertiary'),     # Includes complications (cardio/neuro/disfigurement)
            ss.BoolState('immune'),       # After effective treatment people may acquire temp immunity
            ss.BoolState('ever_exposed'), # Anyone ever exposed - stays true after treatment

            # Congenital syphilis states
            ss.BoolState('congenital'),
            ss.FloatArr('cs_outcome'),

            # Timestep of state changes
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
            ss.FloatArr('ti_normal'),
        )

        return

    @staticmethod
    def make_init_prev_fn(module, sim, uids):
        return sti.make_init_prev_fn(module, sim, uids, active=True)

    @staticmethod
    def make_init_prev_latent_fn(module, sim, uids):
        return sti.make_init_prev_fn(module, sim, uids, active=True, data=module.init_prev_latent_data)

    @property
    def exposed(self):
        """ Default is that exposure equals infection """
        return self.infected

    @property
    def ti_exposed(self):
        """ Alias for ti_infected """
        return self.ti_infected

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

    def init_post(self):
        """ Make initial cases """
        ss.Module.init_post(self) # Avoid super().init_post() since we create infections here
        initial_active_cases = self.pars.init_prev.filter()
        self.set_prognoses(initial_active_cases)
        still_sus = self.susceptible.uids

        # Natural history for initial latent cases
        initial_latent_cases = self.pars.init_latent_prev.filter(still_sus)
        ti_init_cases = self.pars.dist_ti_init_infected.rvs(initial_latent_cases).astype(int)
        self.set_prognoses(initial_latent_cases, ti=ti_init_cases)
        self.set_secondary_prognoses(initial_latent_cases)
        time_to_tertiary = self.pars.time_to_tertiary.rvs(initial_latent_cases)
        self.ti_tertiary[initial_latent_cases] = self.ti_latent[initial_latent_cases] + rr(time_to_tertiary)

        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        results = [
            ss.Result('n_active', dtype=int, label="Number of active cases"),
            ss.Result('pregnant_prevalence', dtype=float, scale=False, label="Pregnant prevalence"),
            ss.Result('detected_pregnant_prevalence', dtype=float, scale=False, label="ANC prevalence"),
            ss.Result('delivery_prevalence', dtype=float, scale=False, label="Delivery prevalence"),
            ss.Result('active_prevalence', dtype=float, scale=False, label="Active prevalence"),
            ss.Result('new_nnds', dtype=int, label="Neonatal deaths"),
            ss.Result('new_stillborns', dtype=int, label="Stillbirths"),
            ss.Result('new_congenital', dtype=int, label="Congenital cases"),
            ss.Result('new_congenital_deaths', dtype=int, label="Congenital deaths"),
            ss.Result('cum_congenital', dtype=int, label="Cumulative congenital cases"),
            ss.Result('cum_congenital_deaths', dtype=int, label="Cumulative congenital deaths"),
            ss.Result('new_deaths', dtype=int, label="Deaths"),

            # Add fetus testing and treatment results, which might be assembled from numerous interventions
            ss.Result('new_fetus_treated_success', dtype=int, label="Fetal treatment success"),
            ss.Result('new_fetus_treated_unnecessary', dtype=int, label="Fetal overtreatment"),
            ss.Result('new_fetus_treated_failure', dtype=int, label="Fetal treatment failure"),
            ss.Result('new_treated_unnecessary_pregnant', dtype=int, label="Overtreatment pregnant"),
        ]

        # Most results are stored by age and sex
        for sk in self.sex_keys.keys():
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' ({sk.upper()})'
            if skk != '':
                results += [
                    ss.Result(f'active_prevalence{skk}', scale=False, label=f"Active prevalence{skl}"),
                ]

            for ab1,ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab2}-{ab2})'
                results += [
                    ss.Result(f'active_prevalence{ask}', scale=False, label=f"Active prevalence{asl}"),
                ]

        # Add FSW and clients to results:
        if self.store_sw:
            results += [
                ss.Result('prevalence_sw', dtype=float, scale=False, label="Prevalence - FSW"),
                ss.Result('new_infections_sw', dtype=int, label="Infections - FSW"),
                ss.Result('new_infections_not_sw', dtype=int, label="Infections - other F"),
                ss.Result('prevalence_client', dtype=float, scale=False, label="Prevalence - clients"),
                ss.Result('new_infections_client', dtype=int, label="Infections - clients"),
                ss.Result('new_infections_not_client', dtype=int, label="Infections - other M"),
            ]

        # Add risk groups to results
        if self.store_risk_groups:
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result('prevalence_risk_group_' + str(risk_group) + '_' + sex, scale=False),
                        ss.Result('new_infections_risk_group_' + str(risk_group) + '_' + sex, dtype=int),
                    ]

        self.define_results(*results)

        return

    def step_state(self):
        """ Updates to states """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Secondary from primary
        secondary_from_primary = self.primary & (self.ti_secondary <= ti)
        if len(secondary_from_primary.uids) > 0:
            self.secondary[secondary_from_primary] = True
            self.primary[secondary_from_primary] = False
            self.set_secondary_prognoses(secondary_from_primary.uids)

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

        # Set rel_trans
        self.rel_trans[self.primary] = self.pars.rel_trans_primary
        self.rel_trans[self.secondary] = self.pars.rel_trans_secondary
        self.rel_trans[self.tertiary] = self.pars.rel_trans_tertiary
        # Latent rel_trans decays with duration of latent infection
        if len(self.latent.uids) > 0:
            self.set_latent_trans()

            # Set people to early
            lu = self.latent.uids
            is_early = lu[(ti - self.ti_latent[lu]) <= self.dur_early[lu]]
            is_late = lu[(ti - self.ti_latent[lu]) > self.dur_early[lu]]
            self.early[is_early] = True
            self.early[is_late] = False
            self.late[is_early] = False
            self.late[is_late] = True

        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        res = self.results
        ppl = self.sim.people

        n_active = res['n_primary'][ti] + res['n_secondary'][ti]
        adults = (ppl.age >= 15) & (ppl.age < 50)
        sexually_active_adults = adults & self.sim.networks.structuredsexual.active(self.sim.people)

        # Overwrite prevalence so we're always storing prevalence of syphilis among sexually active adults
        self.results['prevalence'][ti] = cond_prob(self.infected, sexually_active_adults)
        # self.results['active_prevalence'][ti] = cond_prob(self.active, sexually_active_adults)
        self.results['n_active'][ti] = n_active

        # Pregnant women prevalence, if present
        if 'pregnancy' in self.sim.demographics.keys():
            preg_prev = cond_prob(self.infected, ppl.pregnancy.pregnant)
            self.results['pregnant_prevalence'][ti] = preg_prev
            self.results['detected_pregnant_prevalence'][ti] = preg_prev * self.pars.anc_detection
            deliv_prev = cond_prob(self.infected, ppl.pregnancy.ti_delivery == ti)
            self.results['delivery_prevalence'][ti] = deliv_prev

        # Congenital results
        self.results['new_nnds'][ti]       = np.count_nonzero(self.ti_nnd == ti)
        self.results['new_stillborns'][ti] = np.count_nonzero(self.ti_stillborn == ti)
        self.results['new_congenital'][ti] = np.count_nonzero(self.ti_congenital == ti)
        self.results['new_congenital_deaths'][ti] = self.results['new_nnds'][ti] + self.results['new_stillborns'][ti]
        self.results['new_deaths'][ti] = np.count_nonzero(self.ti_dead == ti)

        # Add FSW and clients to results:
        if self.store_sw:
            fsw = self.sim.networks.structuredsexual.fsw
            clients = self.sim.networks.structuredsexual.client
            self.results['prevalence_sw'][ti] = cond_prob(self.infected, fsw)
            self.results['new_infections_sw'][ti] = np.count_nonzero((self.ti_infected == ti) & fsw)
            self.results['new_infections_not_sw'][ti] = np.count_nonzero((self.ti_infected == ti) & ~fsw)
            self.results['prevalence_client'][ti] = cond_prob(self.infected, clients)
            self.results['new_infections_client'][ti] = np.count_nonzero((self.ti_infected == ti) & clients)
            self.results['new_infections_not_client'][ti] = np.count_nonzero((self.ti_infected == ti) & ~clients)

        # Add risk groups
        if self.store_risk_groups:
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    prev_denom = (self.sim.networks.structuredsexual.risk_group == risk_group) & (self.sim.people[sex])
                    risk_group_new_inf = (self.ti_infected == ti) & (self.sim.networks.structuredsexual.risk_group == risk_group) & (self.sim.people[sex])
                    if risk_group_new_inf.any():
                        self.results['prevalence_risk_group_' + str(risk_group) + '_' + sex][ti] = cond_prob(self.infected, prev_denom)
                        self.results['new_infections_risk_group_' + str(risk_group) + '_' + sex][ti] = np.count_nonzero(risk_group_new_inf)

        # Results by age and sex
        for pkey, pattr in self.sex_keys.items():
            skk = '' if pkey == '' else f'_{pkey}'

            n_act = self.active & ppl[pattr]
            self.results[f'active_prevalence{skk}'][ti] = cond_prob(self.active, adults & ppl[pattr])

            # Compute age results
            age_results = dict(active_prevalence = div(self.agehist(n_act), self.agehist(ppl[pattr])))

            # Store age results
            for akey, ares in age_results.items():
                ai = 0
                for ab1, ab2 in zip(self.age_bins[:-1], self.age_bins[1:]):
                    ask = f'{skk}_{ab1}_{ab2}'
                    self.results[f'{akey}{ask}'][ti] = ares[ai]
                    ai += 1

        return

    def finalize_results(self):
        super().finalize_results()
        self.results['cum_congenital'][:] = np.cumsum(self.results['new_congenital'])
        self.results['cum_congenital_deaths'][:] = np.cumsum(self.results['new_congenital_deaths'])
        return

    def set_latent_trans(self, ti=None):
        if ti is None: ti = self.ti
        dur_latent = ti - self.ti_latent[self.latent]
        hl = self.pars.rel_trans_latent_half_life
        decay_rate = np.log(2) / hl if ~np.isnan(hl) else 0.
        latent_trans = self.pars.rel_trans_latent * np.exp(-decay_rate * dur_latent)
        self.rel_trans[self.latent] = latent_trans
        return

    def set_prognoses(self, uids, source_uids=None, ti=None):
        """
        Set initial prognoses for adults newly infected with syphilis
        """

        if ti is None:
            ti = self.ti
        else:
            # Check that ti is consistent with uids
            if not (sc.isnumber(ti) or len(ti) == len(uids)):
                errormsg = 'ti for set_prognoses must be int or array of length uids'
                raise ValueError(errormsg)

        # Call super method, which records transmissions
        super().set_prognoses(uids, source_uids)

        # Set initial states upon exposure
        self.susceptible[uids] = False
        self.ever_exposed[uids] = True
        self.primary[uids] = True
        self.infected[uids] = True
        self.ti_primary[uids] = ti
        self.ti_infected[uids] = ti

        # Primary to secondary
        dur_primary = self.pars.dur_primary.rvs(uids)
        self.ti_secondary[uids] = self.ti_primary[uids] + rr(dur_primary)
        self.dur_early[uids] = self.pars.dur_early.rvs(uids)

        return

    def set_secondary_prognoses(self, uids):
        """ Set prognoses for people who have just progressed to secondary infection """
        dur_secondary = self.pars.dur_secondary.rvs(uids)
        self.ti_latent[uids] = self.ti_secondary[uids] + rr(dur_secondary)
        return

    def set_latent_prognoses(self, uids):
        # Reactivators
        will_reactivate = self.pars.p_reactivate.rvs(uids)
        reactivate_uids = uids[will_reactivate]
        if len(reactivate_uids) > 0:
            time_to_reactivate = self.pars.time_to_reactivate.rvs(reactivate_uids)
            self.ti_secondary[reactivate_uids] = self.ti_latent[reactivate_uids] + rr(time_to_reactivate)

        # Latent to tertiary
        nonreactivate_uids = uids[~will_reactivate]
        if len(nonreactivate_uids) > 0:
            is_tertiary = self.pars.p_tertiary.rvs(nonreactivate_uids)
            tertiary_uids = nonreactivate_uids[is_tertiary]
            if len(tertiary_uids) > 0:
                time_to_tertiary = self.pars.time_to_tertiary.rvs(tertiary_uids)
                self.ti_tertiary[tertiary_uids] = self.ti_latent[tertiary_uids] + rr(time_to_tertiary)

                # Tertiary to dead
                will_die = self.pars.p_death.rvs(tertiary_uids)
                dead_uids = tertiary_uids[will_die]
                if len(dead_uids) > 0:
                    time_to_death = self.pars.time_to_death.rvs(dead_uids)
                    self.ti_dead[dead_uids] = self.ti_tertiary[dead_uids] + rr(time_to_death)

        return

    def set_congenital(self, target_uids, source_uids=None):
        """
        Natural history of syphilis for congenital infection
        """
        ti = self.ti
        self.susceptible[target_uids] = False
        new_outcomes = {k:0 for k in self.pars.birth_outcome_keys}

        # Determine outcomes
        for state in ['active', 'early', 'late']:

            source_state_inds = getattr(self, state)[source_uids].nonzero()[-1]
            uids = target_uids[source_state_inds]
            source_state_uids = source_uids[source_state_inds]

            if len(uids) > 0:

                # Birth outcomes must be modified to add probability of susceptible birth
                birth_outcomes = self.pars.birth_outcomes[state]
                assigned_outcomes = birth_outcomes.rvs(uids)
                self.cs_outcome[uids] = assigned_outcomes
                timesteps_til_delivery = self.sim.demographics.pregnancy.ti_delivery - self.ti

                # Schedule events
                for oi, outcome in enumerate(self.pars.birth_outcome_keys):
                    m_uids = source_state_uids[assigned_outcomes == oi]
                    o_uids = uids[assigned_outcomes == oi]
                    if len(o_uids) > 0:
                        ti_outcome = f'ti_{outcome}'
                        vals = getattr(self, ti_outcome)
                        vals[o_uids] = ti + timesteps_til_delivery[m_uids]

                        setattr(self, ti_outcome, vals)
                        new_outcomes[outcome] += len(o_uids)

        # Check that the birth outcomes are mutually exclusive
        if not sum(new_outcomes.values()) == len(target_uids):
            raise ValueError('Birth outcomes are not mutually exclusive')

        # Check that the birth outcomes are not greater than the number of congenital cases
        for o1, out1 in enumerate(self.pars.birth_outcome_keys):
            for o2, out2 in enumerate(self.pars.birth_outcome_keys):
                if o1 != o2:
                    val1 = getattr(self, f'ti_{out1}')
                    val2 = getattr(self, f'ti_{out2}')
                    if (val1.notnan & val2.notnan).any() :
                        raise ValueError(f'Birth outcomes {out1} and {out2} are not mutually exclusive')

        return


