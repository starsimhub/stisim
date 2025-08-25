"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI, BaseSTIPars

__all__ = ['HIV', 'HIVPars']


class HIVPars(BaseSTIPars):
    def __init__(self, **kwargs):
        super().__init__()

        # Natural history without treatment
        self.cd4_start = ss.normal(loc=800, scale=50)
        self.cd4_latent = ss.normal(loc=500, scale=50)
        self.dur_acute = ss.lognorm_ex(ss.months(3), ss.months(1))  # Duration of acute HIV infection
        self.dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3))  # Duration of latent, untreated HIV infection
        self.dur_falling = ss.lognorm_ex(ss.years(3), ss.years(1))  # Duration of late-stage HIV when CD4 counts fall
        self.p_hiv_death = ss.bernoulli(p=0)  # Death from HIV-related complications; set by make_p_hiv_death()
        self.include_aids_deaths = True

        # Transmission
        self.beta = 0  # Placeholder, replaced by network-specific betas
        self.beta_m2f = None
        self.rel_beta_f2m = 0.5
        self.beta_m2c = ss.permonth(0.025)  # Approx 0.2 over the course of the pregnany
        self.rel_trans_acute = ss.normal(loc=6, scale=0.5)  # Increase transmissibility during acute HIV infection
        self.rel_trans_falling = ss.normal(loc=8, scale=0.5)  # Increase transmissibility during late HIV infection
        self.eff_condom = 0.9

        # Initialization
        self.init_prev = ss.bernoulli(p=0.05)
        self.rel_init_prev = 1
        self.init_diagnosed = ss.bernoulli(p=0)
        self.dist_ti_init_infected = ss.uniform(low=-10 * 12, high=-5)
        # dist_ti_init_infected=ss.constant(0),  # Experimented with negative values, but safer to use 0

        # Care seeking
        self.care_seeking = ss.normal(loc=1, scale=0.1)  # Distribution of relative care-seeking behavior
        self.maternal_care_scale = 2  # Factor for scaling up care-seeking behavior during pregnancy

        # Treatment effects
        self.art_cd4_growth = 0.1  # Unitless parameter defining how CD4 reconstitutes after starting ART - used in a logistic growth function
        self.art_efficacy = 0.96  # Efficacy of ART
        self.time_to_art_efficacy = ss.months(6)  # Time to reach full ART efficacy (in months) - linear increase in efficacy
        self.art_cd4_pars = dict(cd4_max=1000, cd4_healthy=500)
        self.dur_on_art = ss.lognorm_ex(ss.years(18), ss.years(5))
        self.dur_post_art = ss.normal()  # Note defined in years!
        self.dur_post_art_scale_factor = 0.1

        self.update(kwargs)
        return


class HIV(BaseSTI):

    def __init__(self, pars=None, init_prev_data=None, **kwargs):
        super().__init__()

        # Parameters
        default_pars = HIVPars()
        self.define_pars(**default_pars)
        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data

        # States
        self.define_states(
            # Natural history
            ss.FloatArr('ti_acute'),
            ss.BoolState('acute'),
            ss.FloatArr('ti_latent'),
            ss.BoolState('latent'),
            ss.FloatArr('ti_falling'),
            ss.BoolState('falling'),
            ss.BoolState('post_art'),  # After stopping ART, CD4 falls linearly until death
            ss.FloatArr('ti_zero'),  # Time of zero CD4 count - generally corresponds to AIDS death
            ss.FloatArr('ti_dead'),  # Time of HIV/AIDS death

            # Care and treatment states
            ss.FloatArr('baseline_care_seeking'),
            ss.FloatArr('care_seeking'),
            ss.BoolState('never_art', default=True),
            ss.BoolState('on_art'),
            ss.FloatArr('ti_art'),
            ss.FloatArr('ti_stop_art'),

            # CD4 states
            ss.FloatArr('cd4'),             # Current CD4 count
            ss.FloatArr('cd4_start'),       # Initial CD4 count for each agent before an infection
            ss.FloatArr('cd4_preart'),      # CD4 immediately before initiating ART
            ss.FloatArr('cd4_latent'),      # CD4 count during latent infection
            ss.FloatArr('cd4_nadir'),       # Lowest CD4
            ss.FloatArr('cd4_potential'),   # Potential CD4 count if continually treated
            ss.FloatArr('cd4_postart'),     # CD4 after stopping ART

            # Knowledge of HIV status
            ss.BoolState('diagnosed'),
            ss.FloatArr('ti_diagnosed'),
        )

        return

    @property
    def include_mtct(self): return 'pregnancy' in self.sim.demographics

    def init_results(self):
        """
        Initialize results
        """
        super().init_results()
        results = [
            ss.Result('new_deaths', dtype=int, label='Deaths'),
            ss.Result('cum_deaths', dtype=int, label='Cumulative deaths'),
            ss.Result('new_diagnoses', dtype=int, label='Diagnoses'),
            ss.Result('cum_diagnoses', dtype=int, label='Cumulative diagnoses'),
            ss.Result('new_agents_on_art', dtype=int, label='New treated'),
            ss.Result('prevalence_15_49', dtype=float, label='Prevalence 15-49', scale=False),
            ss.Result('prevalence_sw', dtype=float, label='FSW prevalence', scale=False),
            ss.Result('new_infections_sw', dtype=int, label='New infections - FSW'),
            ss.Result('new_infections_not_sw', dtype=int, label='New infections - Other F'),
            ss.Result('prevalence_client', dtype=float, label='Client prevalence', scale=False),
            ss.Result('new_infections_client', dtype=int, label='New infections - Clients'),
            ss.Result('new_infections_not_client', dtype=int, label='New infections - Other M'),
            ss.Result('p_on_art', dtype=float, label='Proportion on ART', scale=False),
        ]

        if self.include_mtct:
            results += [ss.Result('n_on_art_pregnant', dtype=int)]

        # Add FSW and clients to results:
        if 'structuredsexual' in self.sim.networks.keys():
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result('prevalence_risk_group_' + str(risk_group) + '_' + sex, scale=False),
                        ss.Result('new_infections_risk_group_' + str(risk_group) + '_' + sex, dtype=int),
                    ]

        self.define_results(*results)

        return

    def init_post(self):
        """ Set states """
        ss.Module.init_post(self)  # Skip the disease init_post() since we create infections in a different way
        # Set initial CD4
        self.init_cd4()
        self.init_care_seeking()

        # Make initial cases, some of which may have occured prior to the sim start
        if self.init_prev_data is not None:
            p_init_infection = self.make_init_prev()
            self.pars.init_prev.set(p_init_infection)  # Set the initial prevalence function
        initial_cases = self.pars.init_prev.filter()
        ti_init_cases = self.pars.dist_ti_init_infected.rvs(initial_cases).astype(int)
        self.set_prognoses(initial_cases, ti=ti_init_cases)
        initial_cases_diagnosed = self.pars.init_diagnosed.filter(initial_cases)
        self.diagnosed[initial_cases_diagnosed] = True
        self.ti_diagnosed[initial_cases_diagnosed] = 0
        return

    # CD4 functions
    def acute_decline(self, uids):
        """ Acute decline in CD4 """
        acute_start = self.ti_acute[uids]   # Time of infection
        acute_end = self.ti_latent[uids]    # Time to latent infection
        acute_dur = acute_end - acute_start  # Total time in acute phase, rounded to timesteps
        cd4_start = self.cd4_start[uids]
        cd4_end = self.cd4_latent[uids]
        per_timestep_decline = sc.safedivide(cd4_start-cd4_end, acute_dur)
        cd4 = self.cd4[uids] - per_timestep_decline
        return cd4

    def falling_decline(self, uids):
        """ Decline in CD4 during late-stage infection, when counts are falling """
        falling_start = self.ti_falling[uids]
        falling_end = self.ti_zero[uids]
        falling_dur = falling_end - falling_start
        time_falling = self.ti - self.ti_falling[uids]
        cd4_start = self.cd4_latent[uids]
        cd4_end = 1  # To avoid divide by zero problems
        per_timestep_decline = sc.safedivide(cd4_start-cd4_end, falling_dur)
        cd4 = np.maximum(cd4_end, self.cd4[uids] - per_timestep_decline)
        # cd4 = np.maximum(0, cd4_start - per_timestep_decline*time_falling)
        return cd4

    def post_art_decline(self, uids):
        """
        Decline in CD4 after going off treatment
        This implementation has the possibly-undesirable feature that a person
        who goes on ART for a year and then off again might have a slightly shorter
        lifespan than if they'd never started treatment.
        """
        cd4_end = 1  # To avoid divide by zero problems

        # Death immediately
        zero_now_inds = (self.ti_zero[uids] <= self.ti).nonzero()[-1]
        zero_later_inds = (self.ti_zero[uids] > self.ti).nonzero()[-1]
        cd4 = self.cd4[uids]
        cd4[zero_now_inds] = cd4_end

        # Death later
        zero_later_uids = uids.remove(uids[zero_now_inds])
        ti_zero = self.ti_zero[zero_later_uids]
        ti_stop_art = self.ti_stop_art[zero_later_uids]
        post_art_dur = ti_zero - ti_stop_art
        time_post_art = self.ti - ti_stop_art
        cd4_start = self.cd4_postart[zero_later_uids]
        if post_art_dur.any() <= 0:
            post_art_dur[post_art_dur <= 0] = 1
            error_msg = 'Post-ART duration is negative'
            raise ValueError(error_msg)
        per_timestep_decline = (cd4_start-cd4_end)/post_art_dur
        cd4[zero_later_inds] = np.maximum(cd4_end, cd4_start - per_timestep_decline*time_post_art)
        return cd4

    def cd4_increase(self, uids):
        """
        Increase CD4 counts for people who are receiving treatment.
        Growth curves are calculated to match EMODs CD4 reconstitution equation for people who initiate treatment
        with a CD4 count of 50 (https://docs.idmod.org/projects/emod-hiv/en/latest/hiv-model-healthcare-systems.html)
        However, here we use a logistic growth function and assume that ART CD4 count depends on CD4 at initiation.
        Sources:

            - https://i-base.info/guides/starting/cd4-increase
            - https://www.sciencedirect.com/science/article/pii/S1876034117302022
            - https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-8-20
        """
        # Calculate time on ART and CD4 prior to starting
        ti_art = self.ti_art[uids]
        cd4_preart = self.cd4_preart[uids]
        dur_art = self.ti - ti_art

        # Extract growth parameters
        growth_rate = self.pars.art_cd4_growth
        cd4_total_gain = self.cd4_potential[uids] - self.cd4_preart[uids]
        cd4_now = 2*cd4_total_gain/(1+np.exp(-dur_art*growth_rate))-cd4_total_gain+cd4_preart  # Concave logistic

        return cd4_now

    @property
    def symptomatic(self):
        return self.infectious

    def make_p_hiv_death(self, uids=None):
        cd4_bins = np.array([1000, 500, 350, 200, 50, 0])
        p_hiv_death = ss.peryear(np.array([0.003, 0.003, 0.005, 0.01, 0.05, 0.300])).to_prob(self.dt)
        return p_hiv_death[np.digitize(self.cd4[uids], cd4_bins)]

    @staticmethod
    def _interpolate(vals: list, t):
        vals = sorted(vals, key=lambda x: x[0])  # Make sure values are sorted
        assert len({x[0] for x in vals}) == len(vals)  # Make sure time points are unique
        return np.interp(t, [x[0] for x in vals], [x[1] for x in vals], left=vals[0][1], right=vals[-1][1])

    def init_cd4(self):
        """
        Set CD4 counts
        """
        uids = ss.uids(np.isnan(self.cd4_start.raw).nonzero()[-1])
        self.cd4_start[uids] = self.pars.cd4_start.rvs(uids)
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        return

    def init_care_seeking(self):
        """
        Set care seeking behavior
        """
        uids = ss.uids(self.care_seeking.isnan)
        self.care_seeking[uids] = self.pars.care_seeking.rvs(uids)
        self.baseline_care_seeking[uids] = sc.dcp(self.care_seeking[uids])  # Copy it so pregnancy can modify it
        return

    def step_state(self):
        """
        Carry out autonomous updates at the start of the timestep (prior to transmission)
        """
        ti = self.ti

        # Set initial CD4 counts for new agents:
        self.init_cd4()

        # Handle care seeking behavior. First, initialize, then adjust depending on pregnancy:
        # increase care-seeking for pregnant women and decrease again after postpartum.
        # This makes it much less likely that pregnant women will stop treatment
        self.init_care_seeking()
        if self.include_mtct:
            pregnant = self.sim.demographics.pregnancy.pregnant
            self.care_seeking[pregnant] = self.baseline_care_seeking[pregnant] * self.pars.maternal_care_scale
            self.care_seeking[~pregnant] = self.baseline_care_seeking[~pregnant]

        # Adjust CD4 counts for people receiving treatment - logarithmic increase
        if self.on_art.any():
            art_uids = self.on_art.uids
            self.cd4[art_uids] = self.cd4_increase(art_uids)

        # Adjust CD4 counts for people who have gone off treatment - linear decline
        if (~self.on_art & ~self.never_art).any():
            off_art_uids = (~self.on_art & ~self.never_art).uids
            self.cd4[off_art_uids] = self.post_art_decline(off_art_uids)

        # Update states for people who have never been on ART (ART removes these)
        # Acute & not on ART
        latent = self.acute & (self.ti_latent <= ti)
        falling = self.latent & (self.ti_falling <= ti)
        self.acute[latent] = False
        self.latent[latent] = True
        self.latent[falling] = False
        self.falling[falling] = True

        # Update CD4 counts
        self.cd4[self.acute.uids] = self.acute_decline(self.acute.uids)
        untreated_latent = self.latent
        self.cd4[untreated_latent.uids] = self.cd4_latent[untreated_latent.uids]
        untreated_falling = self.falling
        if untreated_falling.any():
            self.cd4[untreated_falling.uids] = self.falling_decline(untreated_falling.uids)

        # Update CD4 nadir for anyone not on treatment
        untreated = self.infected & ~self.on_art
        self.cd4_nadir[untreated] = np.minimum(self.cd4_nadir[untreated], self.cd4[untreated])

        # Update transmission
        self.update_transmission()

        # Check CD4
        if np.isnan(self.cd4[self.infected]).any():
            errormsg = 'Invalid entry for CD4'
            raise ValueError(errormsg)

        # Update deaths. We capture deaths from AIDS (i.e., when CD4 count drops to ~0) as well as deaths from
        # serious HIV-related illnesses, which can occur throughout HIV.
        off_art = (self.infected & ~self.on_art).uids
        p_death = self.make_p_hiv_death(uids=off_art)
        self.pars.p_hiv_death.set(0)
        self.pars.p_hiv_death.set(p_death)  # Set the death probability function
        hiv_deaths = self.pars.p_hiv_death.filter(off_art)
        if len(hiv_deaths):
            self.ti_dead[hiv_deaths] = ti
            self.sim.people.request_death(hiv_deaths)
        if self.pars.include_aids_deaths:
            aids_deaths = (self.ti_zero <= ti).uids
            if len(aids_deaths):
                self.ti_dead[aids_deaths] = ti
                self.sim.people.request_death(aids_deaths)
        return

    def step_die(self, uids):
        """ Clear all states for dead agents """

        # Reset boolean states
        self.susceptible[uids] = False
        self.infected[uids] = False
        self.acute[uids] = False
        self.latent[uids] = False
        self.falling[uids] = False
        self.post_art[uids] = False
        self.never_art[uids] = False
        self.on_art[uids] = False
        self.diagnosed[uids] = False

        # Clear time states except for ti_dead
        self.ti_infected[uids] = np.nan
        self.ti_acute[uids] = np.nan
        self.ti_latent[uids] = np.nan
        self.ti_falling[uids] = np.nan
        self.ti_zero[uids] = np.nan
        self.ti_art[uids] = np.nan
        self.ti_stop_art[uids] = np.nan
        self.ti_diagnosed[uids] = np.nan

        # Clear CD4 states
        self.cd4[uids] = np.nan
        self.cd4_start[uids] = np.nan
        self.cd4_preart[uids] = np.nan
        self.cd4_latent[uids] = np.nan
        self.cd4_nadir[uids] = np.nan
        self.cd4_potential[uids] = np.nan
        self.cd4_postart[uids] = np.nan

        # Clear all other states
        self.rel_sus[uids] = np.nan
        self.rel_trans[uids] = np.nan
        self.baseline_care_seeking[uids] = np.nan
        self.care_seeking[uids] = np.nan

        return

    def update_transmission(self):
        """
        Update rel_trans and rel_sus for all agents. These are reset on each timestep then adjusted depending on states.
        Adjustments are made throughout different modules:

           - rel_trans for acute and late-stage untreated infection are adjusted below
           - rel_trans for all people on treatment (including pregnant women) below
           - rel_sus for unborn babies of pregnant WLHIV receiving treatment is adjusted in the ART intervention
        """
        ti = self.ti

        # Reset susceptibility and infectiousness
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Update rel_trans to account for acute and late-stage infection
        self.rel_trans[self.acute] *= self.pars.rel_trans_acute.rvs(self.acute.uids)
        aids = self.cd4 < 200
        self.rel_trans[aids] *= self.pars.rel_trans_falling.rvs(aids.uids)

        # Update transmission for agents on ART
        # When agents start ART, determine the reduction of transmission (linearly decreasing over 6 months)
        if self.on_art.any():
            full_eff = self.pars.art_efficacy
            time_to_full_eff = self.pars.time_to_art_efficacy
            art_uids = self.on_art.uids
            timesteps_on_art = ti - self.ti_art[art_uids]
            new_on_art = timesteps_on_art < time_to_full_eff/self.dt
            efficacy_to_date = np.full_like(art_uids, fill_value=full_eff, dtype=float)
            efficacy_to_date[new_on_art] = timesteps_on_art[new_on_art]*full_eff/time_to_full_eff.value
            self.rel_trans[art_uids] *= 1 - efficacy_to_date

        return

    def update_results(self):
        """
        Update results at each time step
        """
        super().update_results()
        ti = self.ti

        # Recalculate prevalence so it's for the whole population - the STI module calculates it for adults
        self.results['prevalence'][ti] = sum(self.infected) / len(self.infected)

        self.results['new_deaths'][ti] = np.count_nonzero(self.ti_dead == ti)
        self.results['cum_deaths'][ti] = np.sum(self.results['new_deaths'][:ti + 1])
        self.results['new_diagnoses'][ti] = np.count_nonzero(self.ti_diagnosed == ti)
        self.results['cum_diagnoses'][ti] = np.sum(self.results['new_diagnoses'][:ti + 1])
        self.results['new_agents_on_art'][ti] = np.count_nonzero(self.ti_art == ti)
        if self.include_mtct:
            self.results['n_on_art_pregnant'][ti] = np.count_nonzero(self.on_art & self.sim.people.pregnancy.pregnant)
        self.results['p_on_art'][ti] = sc.safedivide(self.results['n_on_art'][ti], self.results['n_infected'][ti])

        # Subset by age group:
        infected_15_19 = self.infected[(self.sim.people.age >= 15) & (self.sim.people.age < 50)]
        self.results['prevalence_15_49'][ti] = sum(infected_15_19) / len(infected_15_19)

        # Subset by FSW and client:
        if 'structuredsexual' in self.sim.networks.keys():
            fsw_infected = self.infected[self.sim.networks.structuredsexual.fsw]
            client_infected = self.infected[self.sim.networks.structuredsexual.client]

            # Add FSW and clients to results:
            if len(fsw_infected) > 0:
                self.results['prevalence_sw'][ti] = sum(fsw_infected) / len(fsw_infected)
                self.results['new_infections_sw'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.fsw).uids)
                self.results['new_infections_not_sw'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.fsw).uids)
            if len(client_infected) > 0:
                self.results['prevalence_client'][ti] = sum(client_infected) / len(client_infected)
                self.results['new_infections_client'][ti] = len(((self.ti_infected == ti) & self.sim.networks.structuredsexual.client).uids)
                self.results['new_infections_not_client'][ti] = len(((self.ti_infected == ti) & ~self.sim.networks.structuredsexual.client).uids)

        return

    def set_prognoses(self, uids, sources=None, ti=None):
        """
        Set prognoses upon infection
        """
        if ti is None:
            ti = self.ti
        else:
            # Check that ti is consistent with uids
            if not (sc.isnumber(ti) or len(ti) == len(uids)):
                errormsg = 'ti for set_prognoses must be int or array of length uids'
                raise ValueError(errormsg)

        self.susceptible[uids] = False
        self.infected[uids] = True
        self.acute[uids] = True

        self.ti_infected[uids] = ti
        self.ti_acute[uids] = ti
        self.cd4[uids] = self.cd4_start[uids]

        # Set timing and CD4 count of latent infection
        dur_acute = self.pars.dur_acute.rvs(uids)
        self.ti_latent[uids] = self.ti_acute[uids] + dur_acute.astype(int)
        self.cd4_latent[uids] = self.pars.cd4_latent.rvs(uids)

        # Set time of onset of late-stage CD4 decline
        dur_latent = self.pars.dur_latent.rvs(uids)
        self.ti_falling[uids] = self.ti_latent[uids] + dur_latent.astype(int)
        dur_falling = self.pars.dur_falling.rvs(uids)
        self.ti_zero[uids] = self.ti_falling[uids] + dur_falling.astype(int)

        return

    def set_congenital(self, uids, sources):
        self.cd4_start[uids] = sc.dcp(self.cd4_start[sources])
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        self.set_prognoses(uids, sources)
        return

    # Treatment-related changes
    def start_art(self, uids):
        """
        Check who is ready to start ART treatment and put them on ART
        """
        ti = self.ti

        self.on_art[uids] = True
        newly_treated = uids[self.never_art[uids]]
        self.never_art[newly_treated] = False
        self.ti_art[uids] = ti
        self.cd4_preart[uids] = self.cd4[uids]

        # Determine when agents goes off ART
        dur_on_art = self.pars.dur_on_art.rvs(uids)
        self.ti_stop_art[uids] = ti + dur_on_art.astype(int)

        # ART nullifies all states and all future dates in the natural history
        self.acute[uids] = False
        self.latent[uids] = False
        self.falling[uids] = False
        future_latent = uids[self.ti_latent[uids] > ti]
        self.ti_latent[future_latent] = np.nan
        future_falling = uids[self.ti_falling[uids] > ti]
        self.ti_falling[future_falling] = np.nan
        future_zero = uids[self.ti_zero[uids] > ti]  # NB, if they are scheduled to die on this time step, they will
        self.ti_zero[future_zero] = np.nan

        # Set CD4 potential for anyone new to treatment - retreated people have the same potential
        # Extract growth parameters
        if len(newly_treated) > 0:
            cd4_max = self.pars.art_cd4_pars['cd4_max']
            cd4_healthy = self.pars.art_cd4_pars['cd4_healthy']
            cd4_preart = self.cd4_preart[newly_treated]

            # Calculate potential CD4 increase - assuming that growth follows the concave part of a logistic function
            # and that the total gain depends on the CD4 count at initiation
            if (cd4_preart == 0).any():
                raise ValueError('CD4 count is zero at initiation of ART')
            cd4_scale_factor = (cd4_max-cd4_preart)/cd4_healthy*np.log(cd4_max/cd4_preart)
            cd4_total_gain = cd4_preart*cd4_scale_factor
            self.cd4_potential[newly_treated] = self.cd4_preart[newly_treated] + cd4_total_gain

        return

    def stop_art(self, uids=None):
        """
        Check who is stopping ART treatment and put them off ART
        """
        ti = self.ti

        # Remove agents from ART
        if uids is None: uids = self.on_art & (self.ti_stop_art <= ti)
        self.on_art[uids] = False
        self.ti_stop_art[uids] = ti
        self.cd4_postart[uids] = sc.dcp(self.cd4[uids])

        # Set decline
        dur_mean = np.log(self.cd4_preart[uids])*self.cd4[uids]/self.cd4_potential[uids]
        dur_scale = dur_mean * self.pars.dur_post_art_scale_factor
        dur_mean = ss.years(dur_mean)
        dur_scale = ss.years(np.maximum(dur_scale, 1e-3))  # Ensure it's not zero
        self.pars.dur_post_art.set(loc=dur_mean, scale=dur_scale)
        dur_post_art = self.pars.dur_post_art.rvs(uids)
        if np.isnan(dur_post_art).any():
            errormsg = 'Invalid entry for post-ART duration'
            raise ValueError(errormsg)
        self.ti_zero[uids] = ti + dur_post_art.astype(int)

        return

