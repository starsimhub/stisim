"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
from stisim.diseases.sti import STIResults, init_prev_by_group, sexual_beta, hetero_networks, default_age_bins, default_sex_keys

__all__ = ['HIV']

# The parameters of sti.ART, which HIV takes and passes on to it (see sti.ART.init_pre()), as STIsim 3 gave them to HIV
art_pars = ['rel_art_mortality_effective', 'rel_art_mortality_unsupp_m', 'rel_art_mortality_unsupp_f', 'art_death_age', 'art_death_dur',
            'art_cd4_growth', 'art_cd4_pars', 'care_seeking', 'maternal_care_scale', 'effective_art_efficacy', 'nonsupp_art_efficacy',
            'rate_art_dropout', 'time_to_art_efficacy', 'p_effective_art', 'dur_on_art', 'rel_dur_on_art', 'dur_on_art_trend',
            'dur_post_art', 'dur_post_art_scale_factor']

# The states and methods of sti.ART that can be read through HIV (e.g. hiv.on_effective_art), as in STIsim 3
art_attrs = ['on_effective_art', 'on_nonsuppressive_art', 'art_discontinued', 'art_naive', 'ti_art', 'ti_stop_art', 'cd4_preart',
             'cd4_potential', 'cd4_postart', 'dur_post_art', 'care_seeking', 'baseline_care_seeking', 'start_art', 'stop_art',
             'set_suppression', 'art_factor', 'on_art_death_rate']


class HIV(ss.Infection, STIResults):
    """
    HIV: its natural history, CD4 count, transmission, diagnosis, PrEP, and circumcision

    People infected are in acute infection for `dur_acute`, then latent for `dur_latent`, then have a falling
    CD4 count for `dur_falling`, after which they die of AIDS (with probability `p_aids_death`); these times
    are drawn at infection. Throughout, they die of HIV-related illness at a rate by CD4 count (`cd4_death_rates`).
    Infection is by sexual transmission (`beta_m2f`), and from mother to child in pregnancy (`beta_m2c`) and
    breastfeeding (`beta_breastfeed`). The CD4 count falls in acute and late infection (by the share of the stage
    elapsed). Infectiousness is by stage, and susceptibility by PrEP, circumcision, and age.

    Treatment is a module of its own (`sti.ART`), which registers itself with `add_treatment()`. HIV asks each
    treatment, at set points in its step, for the CD4 counts of the people it treats, and for multipliers of their
    death rates and infectiousness; a treatment takes people out of the natural history when it starts them (see
    `sti.ART`). `on_art` and the other states of treatment (`on_effective_art`, `ti_art`, ...) are read from it.
    Testing refers people to treatment with `refer()`.

    States: `susceptible`, `infected`, `acute`, `latent`, and `falling`; `on_art` and `aids` are conditions.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'hiv'``)
        init_prev_data (DataFrame): initial prevalence by risk group, sex, and sex work, in place of `init_prev`
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care
        kwargs (dict): parameter values (e.g. ``beta_m2f=0.02``), `init_prev`, or `label`
    """
    hetero_networks = hetero_networks # The sexual networks between men and women, by name, which a subclass can add to

    def __init__(self, pars=None, name='hiv', init_prev_data=None, age_bins=default_age_bins,
                 sex_keys=default_sex_keys, include_care=True, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            dt = 'month',
            init_prev = 0.05,

            # Transmission: a probability per act on the sexual networks, and rates per month from mother to child
            beta = 0.0, # On each network, from the ones below (see validate_beta())
            beta_m2f = 0.05, # From a man to a woman
            rel_beta_f2m = 0.5, # From a woman to a man, relative to from a man to a woman
            beta_m2m = 0.0, # Between men
            beta_m2c = 0.025, # From mother to child in pregnancy, per month (about 20% over a pregnancy)
            beta_breastfeed = 0.005, # From mother to child by breastfeeding, per month (about 14% over 12 months without ART)
            eff_condom = 0.9, # Read by the sexual networks (see sti.BaseNetwork.net_beta())

            # Natural history
            dur_acute = ss.lognorm_ex(ss.months(1.7), ss.months(1)), # Bellan 2015: central estimate 1.7 months
            dur_latent = ss.lognorm_ex(ss.years(10), ss.years(3)), # Untreated
            dur_falling = ss.lognorm_ex(ss.years(3), ss.years(1)), # Late-stage infection, as the CD4 count falls to zero
            p_aids_death = 1.0, # Of dying of AIDS when the CD4 count reaches zero

            # Death from HIV-related illness, by CD4 count, off ART
            cd4_death_bins = np.array([1000, 500, 350, 200, 50, 0]), # Off-ART CD4-stratified mortality (descending)
            cd4_death_rates = np.array([0.003, 0.003, 0.005, 0.01, 0.05, 0.30]), # Annual rate per bin
            rel_death = 1.0, # Scales all HIV death rates (off- and on-ART)
            rel_death_f = 0.74, # Additional mortality multiplier for females, on- and off-ART alike

            # CD4 count
            cd4_start = ss.normal(loc=800, scale=50),
            cd4_latent = ss.normal(loc=500, scale=50),

            # Infectiousness and susceptibility
            rel_trans_acute = ss.normal(loc=5.3, scale=0.5), # Bellan 2015: RH=5.3, EHM≈7.3
            rel_trans_falling = ss.normal(loc=8, scale=0.5), # Increase in late infection (a CD4 count below 200)
            rel_sus_age = None, # Age/sex-stratified rel_sus multipliers: (age_lo, age_hi, sex, multiplier) tuples; sex is 'f', 'm', or None for both. E.g. [(15, 25, 'f', 1.7), (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)]
            eff_circ = 0.6, # Program (medical) VMMC: reduction in HIV acquisition susceptibility
            eff_circ_traditional = 0.6, # Traditional (non-program) circumcision: reduction in HIV acquisition susceptibility

            # Initial conditions
            rel_init_prev = 1, # Scales init_prev_data
            init_diagnosed = ss.bernoulli(p=0),
            dist_ti_init_infected = ss.uniform(low=-10*12, high=-5), # When the infections at the start began, in steps

            # Treatment: sti.ART's parameters, which it takes from here if they are given to HIV (as in STIsim 3)
            **{key: None for key in art_pars},
        )
        self.update_pars(pars, **kwargs)
        self.init_prev_data = init_prev_data
        if init_prev_data is not None:
            self.pars.init_prev = ss.bernoulli(p=self.init_prev_by_group)

        self.define_states(
            # Natural history and treatment
            ss.BoolState('acute'),
            ss.BoolState('latent'),
            ss.BoolState('falling'),
            ss.TimeArr('ti_latent'), # When each stage starts, drawn at infection
            ss.TimeArr('ti_falling'),
            ss.TimeArr('ti_zero'), # When the CD4 count reaches zero, at the end of late infection or after stopping treatment: an AIDS death
            ss.FloatArr('ti_dead'), # Time of death of HIV

            # CD4 count
            ss.FloatArr('cd4'), # Current CD4 count
            ss.FloatArr('cd4_start'), # Initial CD4 count for each agent before an infection
            ss.FloatArr('cd4_latent'), # CD4 count during latent infection
            ss.FloatArr('cd4_nadir'), # Lowest CD4

            # Knowledge of HIV status
            ss.BoolState('diagnosed'),
            ss.FloatArr('ti_diagnosed'),

            # PrEP. Efficacy/duration/adherence are specified directly on whichever Prep intervention enrolls someone:
            # prep_eff is the realized efficacy for this person's current course (eff*adherence, computed at enrollment),
            # and prep_source identifies which Prep instance enrolled them (for its own coverage-target accounting)
            ss.BoolState('prep_naive', default=True), # Never started PrEP, from any source
            ss.BoolState('on_prep'),
            ss.BoolState('prep_discontinued'), # Currently off after a course ended
            ss.FloatArr('prep_eff'), # Realized rel_sus multiplier reduction for the current course
            ss.FloatArr('prep_source'), # Id of the enrolling Prep instance; NaN=not on PrEP
            ss.FloatArr('ti_prep_start'),
            ss.TimeArr('ti_prep_stop'), # Scheduled course expiry (renewal/lapse), like ti_stop_art

            # Circumcision (VMMC)
            ss.BoolState('circumcised'),
            ss.BoolState('circ_traditional', default=False), # True=traditional (non-program), False=program (medical) VMMC
            ss.FloatArr('ti_circumcised'),
        )
        self.u_death = ss.random(name='u_death') # One uniform a step for each person infected, for death
        self.treatment_names = [] # The modules that treat HIV (e.g. sti.ART), which register themselves (see add_treatment()), by name, so that their random numbers stay their own
        self.set_strata([15, 50], age_bins, sex_keys)
        self.include_care = include_care
        return

    # Conditions over the states
    @property
    def symptomatic(self):
        return self.infectious

    @property
    def treatable(self):
        return self.infected

    @property
    def on_art(self):
        """ On treatment (e.g. ART, suppressed or not), by any treatment """
        out = self.infected & ~self.infected
        for tx in self.treatments:
            out = out | tx.on_treatment
        return out

    # Treatment: the interface for the modules that treat HIV (e.g. sti.ART)
    def add_treatment(self, treatment):
        """
        Register a module that treats HIV. It provides `on_treatment` (a state), and is called at set points in
        HIV's step: `update_cd4()` sets the CD4 counts of those it treats (or has treated), `step_treatment()`
        runs before deaths (e.g. stopping early), `death_multiplier(uids)` and `trans_multiplier(uids)` scale the
        death rates and infectiousness of those on it, `refer(uids, ti)` schedules a start after a diagnosis, and
        `step_die(uids)` clears those who die
        """
        if treatment.name not in self.treatment_names:
            self.treatment_names.append(treatment.name)
        return

    @property
    def treatments(self):
        """ The modules that treat HIV """
        return [self.sim.interventions[name] for name in self.treatment_names] if self.treatment_names else []

    @property
    def art(self):
        """ The treatment that holds the states of ART (`sti.ART`), if any """
        return next((tx for tx in self.treatments if hasattr(tx, 'on_effective_art')), None)

    def refer(self, uids, ti=None):
        """ Refer people diagnosed to treatment, to start at `ti` (by default now) """
        for tx in self.treatments:
            tx.refer(uids, self.ti if ti is None else ti)
        return

    def death_multiplier(self, uids):
        """ The treatments' multiplier of the death rates of some infected agents (1 for those not on treatment) """
        mult = np.ones(len(uids))
        for tx in self.treatments:
            mult *= tx.death_multiplier(uids)
        return mult

    @property
    def aids(self):
        return self.cd4 < 200

    @property
    def include_mtct(self):
        return 'pregnancy' in self.sim.demographics

    def init_prev_by_group(self, sim, uids):
        """ The initial prevalence of each agent, from init_prev_data """
        return init_prev_by_group(sim, uids, self.init_prev_data, self.pars.rel_init_prev)

    def validate_beta(self):
        """ The betas of the sexual networks and the maternal ones, from the parameters (so a fitted beta_m2f is read as it is now) """
        p = self.pars
        sexual = sexual_beta(p.beta_m2f, p.beta_m2f*p.rel_beta_f2m, p.beta_m2m, hetero=self.hetero_networks)
        self.pars.beta = sexual | dict(maternalnet=ss.permonth(p.beta_m2c), breastfeedingnet=ss.permonth(p.beta_breastfeed))
        return super().validate_beta()

    def plot(self):
        """ Plot key HIV results """
        import matplotlib.pyplot as plt
        res = self.results
        keys = ['new_infections', 'new_deaths', 'prevalence', 'prevalence_15_49', 'new_diagnoses', 'p_on_art']
        keys = [k for k in keys if k in res] # Only plot keys that exist
        with sc.options.with_style('fancy'):
            fig, axs = sc.getrowscols(len(keys), make=True)
            for ax, k in zip(axs.flatten(), keys):
                ax.plot(res.timevec, res[k], lw=2, alpha=0.8)
                ax.set_title(res[k].label)
                ax.set_xlabel('Year')
            sc.figlayout()
        return ss.return_fig(fig)

    # Death
    def off_art_death_rate(self, uids):
        """ The annual rate of death from HIV-related illness of each agent off ART, by CD4 count, sex, and rel_death """
        rate = self.pars.cd4_death_rates[np.digitize(self.cd4[uids], self.pars.cd4_death_bins)] * self.pars.rel_death
        rate[~self.sim.people.male[uids]] *= self.pars.rel_death_f
        return rate

    def death_hazard(self, uids):
        """ The annual hazard of death from HIV-related illness of each agent, by CD4 count, lower on ART (not the AIDS death when the CD4 count reaches zero, at a set time), e.g. for MIGHTI's competing risks of death """
        uids = ss.uids(uids)
        h = np.zeros(len(uids))
        infected = uids[self.infected[uids]]
        h[self.infected[uids]] = self.off_art_death_rate(infected)*self.death_multiplier(infected)
        return h

    def init_post(self):
        """ CD4 counts and care seeking; and the infections at the start, which began before it (`dist_ti_init_infected`) """
        self.init_cd4()
        ss.Disease.init_post(self) # Not ss.Infection's, which seeds the infections at the start of the run
        ppl = self.sim.people
        seeds = self.pars.init_prev.filter(ppl.auids[ppl.age.raw[ppl.auids] >= 0]) # Not the fetuses of pregnancies under way
        self.set_prognoses(seeds, ti=self.pars.dist_ti_init_infected.rvs(seeds).astype(int))
        diagnosed = self.pars.init_diagnosed.filter(seeds)
        self.diagnosed[diagnosed] = True
        self.ti_diagnosed[diagnosed] = 0
        self.refer(diagnosed, ti=0) # To start ART at once
        unused = [key for key in art_pars if self.pars[key] is not None] # Taken by sti.ART when there is one (see its init_pre())
        if unused:
            errormsg = f'ART parameters were given to HIV ({sc.strjoin(unused)}), but there is no ART in the sim to use them: add sti.ART()'
            raise ValueError(errormsg)
        return

    # CD4 count and care seeking
    def init_cd4(self):
        """ Set the CD4 counts of new agents """
        uids = ss.uids(np.isnan(self.cd4_start.raw).nonzero()[-1])
        self.cd4_start[uids] = self.pars.cd4_start.rvs(uids)
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        return

    def update_cd4(self):
        """ CD4 counts: falling from the start of infection to the latent level in acute infection, and from it to 1 in late infection; and those of the people treatments have taken over, from them """
        ti = self.ti
        share = lambda start, end: np.clip(np.divide(ti - start, end - start, out=np.ones(len(start)), where=end > start), 0, 1) # Of a stage elapsed
        acute = self.acute.uids
        cd4_start, cd4_latent = self.cd4_start[acute], self.cd4_latent[acute]
        self.cd4[acute] = cd4_start - (cd4_start - cd4_latent)*share(self.ti_infected[acute], self.ti_latent[acute])
        latent = self.latent.uids
        self.cd4[latent] = self.cd4_latent[latent]
        falling = self.falling.uids
        self.cd4[falling] = np.maximum(1, self.cd4_latent[falling] - (self.cd4_latent[falling] - 1)*share(self.ti_falling[falling], self.ti_zero[falling]))
        for tx in self.treatments:
            tx.update_cd4()
        untreated = (self.infected & ~self.on_art).uids
        self.cd4_nadir[untreated] = np.minimum(self.cd4_nadir[untreated], self.cd4[untreated])
        return

    def update_transmission(self):
        """
        Each agent's susceptibility (age, PrEP, circumcision) and infectiousness (acute and late infection, and ART, which
        reaches its full efficacy by a linear rise over time_to_art_efficacy), reset each step; connectors change them after this
        """
        p = self.pars
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        # Age/sex-stratified susceptibility, multiplied in alongside the other contributors, so values are ceteris paribus
        if p.rel_sus_age is not None:
            ppl = self.sim.people
            for age_lo, age_hi, sex, mult in p.rel_sus_age:
                sex_mask = ppl.female if sex == 'f' else ppl.male if sex == 'm' else ppl.alive
                self.rel_sus[sex_mask & (ppl.age >= age_lo) & (ppl.age < age_hi)] *= mult

        # PrEP and circumcision reduce susceptibility
        if self.on_prep.any():
            self.rel_sus[self.on_prep] *= 1 - self.prep_eff[self.on_prep]
        if self.circumcised.any():
            self.rel_sus[self.circumcised & ~self.circ_traditional] *= 1 - p.eff_circ
            self.rel_sus[self.circumcised & self.circ_traditional] *= 1 - p.eff_circ_traditional

        # Infectiousness in acute and late infection, and on ART
        self.rel_trans[self.acute] *= p.rel_trans_acute.rvs(self.acute.uids)
        self.rel_trans[self.aids] *= p.rel_trans_falling.rvs(self.aids.uids)
        for tx in self.treatments:
            treated = tx.on_treatment.uids
            if len(treated):
                self.rel_trans[treated] *= tx.trans_multiplier(treated)
        return

    def infectiousness(self, uids=None):
        """ How infectious each agent is, all told, relative to untreated latent infection (their rel_trans) """
        return self.rel_trans[self.sim.people.auids if uids is None else ss.uids(uids)]

    def step_state(self):
        """ CD4 counts; the natural history; treatment's own steps; death; then each agent's susceptibility and infectiousness """
        ti = self.ti
        self.init_cd4()

        # The natural history, untreated: treatment takes people out of it (see sti.ART.start_art())
        latent = (self.acute & self.ti_latent.due()).uids
        self.acute[latent] = False
        self.latent[latent] = True
        falling = (self.latent & self.ti_falling.due()).uids
        self.latent[falling] = False
        self.falling[falling] = True
        self.update_cd4()

        # Treatment's own steps (e.g. stopping early)
        for tx in self.treatments:
            tx.step_treatment()

        # Death: of AIDS once the CD4 count reaches zero, and at any time at a rate by CD4 count, scaled by treatment
        infected = self.infected.uids
        on_art = self.on_art[infected]
        p_death = ss.peryear(self.off_art_death_rate(infected)*self.death_multiplier(infected)).to_prob(self.dt)
        zero = self.ti_zero.due().raw[infected] & ~on_art
        p_death[zero] = 1 - (1 - p_death[zero])*(1 - self.pars.p_aids_death)
        died = infected[self.u_death.rvs(infected) < p_death]
        self.ti_zero[infected[zero]] = np.nan # Those who survive it live on, at the rates by CD4 count
        if len(died):
            self.ti_dead[died] = ti
            self.sim.people.request_death(died)

        self.update_transmission()
        return

    def step(self):
        """ Transmission (a child infected by its mother starting with her CD4 count), and the new infections by route """
        new, sources, networks = self.transmit()
        nets = list(self.sim.networks.values())
        prenatal = np.isin(networks, [i for i,net in enumerate(nets) if isinstance(net, ss.PrenatalNet)])
        postnatal = np.isin(networks, [i for i,net in enumerate(nets) if isinstance(net, ss.PostnatalNet)])
        mtct = prenatal | postnatal
        if len(new):
            self.cd4_start[new[mtct]] = self.cd4_start[sources[mtct]]
            self.cd4_nadir[new[mtct]] = self.cd4_start[new[mtct]]
            self.set_prognoses(new, sources)
            if self.infection_log:
                self.infection_log.add_data(new, network=networks)
        ti = self.ti
        res = self.results
        res.new_infections_sex[ti] = np.count_nonzero(~mtct)
        res.new_infections_prenatal[ti] = np.count_nonzero(prenatal)
        res.new_infections_postnatal[ti] = np.count_nonzero(postnatal)
        res.new_infections_mtct[ti] = np.count_nonzero(mtct)
        self.sexual_sources = sources[~mtct]
        return new, sources, networks

    def set_prognoses(self, uids, sources=None, ti=None):
        """ Infect some people (at `ti`, by default now; one for each for the infections at the start), and draw their natural history: when each stage starts """
        super().set_prognoses(uids, sources)
        start = (self.ti + 1 if self.sim.initialized else self.ti) if ti is None else ti # A course runs from when the infection happens: at the end of the step, in transmission (as a declared disease dates it), or when the seeds at the start were infected
        ti = self.ti if ti is None else ti
        self.infected[uids] = True
        self.acute[uids] = True
        self.ti_infected[uids] = ti
        self.cd4[uids] = self.cd4_start[uids]
        self.cd4_latent[uids] = self.pars.cd4_latent.rvs(uids)
        p = self.pars
        self.ti_latent[uids] = start + p.dur_acute.rvs(uids) # Due at the nearest step (see ss.TimeArr.due())
        self.ti_falling[uids] = self.ti_latent[uids] + p.dur_latent.rvs(uids)
        self.ti_zero[uids] = self.ti_falling[uids] + p.dur_falling.rvs(uids)
        return

    def step_die(self, uids):
        """ The people who died leave the states of HIV """
        for state in [self.infected, self.acute, self.latent, self.falling, self.on_prep, self.susceptible]:
            state[uids] = False
        for tx in self.treatments:
            tx.step_die(uids)
        return

    def new_infected(self):
        """ Who was infected in this step (not those at the start, whose infections began before it) """
        return self.ti_infected == self.ti


    def start_prep(self, uids, eff, dur, source_id, adh=1.0):
        """
        Start a PrEP course for uids not already on PrEP

        Args:
            uids: agents to start
            eff: base efficacy (0-1) of this course when fully adherent
            dur: course duration (an ss.dur/Time value) before renewal is needed
            source_id: numeric id of the enrolling Prep instance (for that instance's own coverage-target accounting)
            adh: adherence level (0-1), multiplied into eff to get the realized efficacy (default 1.0 = fully adherent)

        Returns:
            The subset of uids actually started (already-on-PrEP agents excluded).
        """
        uids = uids[~self.on_prep[uids]]
        if len(uids) == 0:
            return uids
        ti = self.ti
        self.on_prep[uids] = True
        self.prep_discontinued[uids] = False # Re-starting PrEP clears the discontinued flag
        newly_started = uids[self.prep_naive[uids]]
        self.prep_naive[newly_started] = False
        self.prep_eff[uids] = eff * adh
        self.prep_source[uids] = source_id
        self.ti_prep_start[uids] = ti
        self.ti_prep_stop[uids] = ti + int(dur / self.dt)
        return uids

    def stop_prep(self, uids=None):
        """ Stop PrEP for uids, or (if None) everyone whose current course has expired """
        if uids is None:
            uids = (self.on_prep & self.ti_prep_stop.due()).uids
        self.on_prep[uids] = False
        self.prep_discontinued[uids] = True
        self.prep_eff[uids] = np.nan
        self.prep_source[uids] = np.nan
        self.ti_prep_start[uids] = np.nan
        self.ti_prep_stop[uids] = np.nan
        return uids

    def circumcise(self, uids, traditional=False):
        """
        Mark agents as circumcised. Idempotent: already-circumcised agents in ``uids`` keep their original
        ``circ_traditional``/``ti_circumcised``, so a traditional-MC pass can't overwrite someone circumcised by a program.

        Args:
            uids: agents to circumcise
            traditional: False=program (medical) VMMC (default), True=traditional (non-program)
        """
        uids = uids[~self.circumcised[uids]]
        self.circumcised[uids] = True
        self.circ_traditional[uids] = traditional
        self.ti_circumcised[uids] = self.ti
        return

    # Results
    def init_results(self):
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        results = [
            ss.Result('on_art', dtype=int, kind='stock', label='On ART'),
            ss.Result('new_deaths', dtype=int, label='Deaths'),
            ss.Result('cum_deaths', dtype=int, label='Cumulative deaths', auto_plot=False),
            ss.Result('new_diagnoses', dtype=int, label='Diagnoses'),
            ss.Result('cum_diagnoses', dtype=int, label='Cumulative diagnoses', auto_plot=False),
            ss.Result('new_agents_on_art', dtype=int, label='New treated', auto_plot=False),
            ss.Result('prevalence_15_49', dtype=float, label='Prevalence 15-49', scale=False),
            ss.Result('prevalence_sw', dtype=float, label='FSW prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_sw', dtype=int, label='New infections - FSW', auto_plot=False),
            ss.Result('new_infections_not_sw', dtype=int, label='New infections - Other F', auto_plot=False),
            ss.Result('prevalence_client', dtype=float, label='Client prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_client', dtype=int, label='New infections - Clients', auto_plot=False),
            ss.Result('new_infections_not_client', dtype=int, label='New infections - Other M', auto_plot=False),
            ss.Result('p_on_art', dtype=float, label='Proportion on ART', scale=False),
        ]
        if self.include_mtct:
            results += [
                ss.Result('n_on_art_pregnant', dtype=int, auto_plot=False),
                ss.Result('p_diagnosed_pregnant', dtype=float, label='Proportion of HIV+ pregnant women diagnosed', scale=False, auto_plot=False),
            ]
        if 'structuredsexual' in self.sim.networks.keys():
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result('prevalence_risk_group_' + str(risk_group) + '_' + sex, scale=False, auto_plot=False),
                        ss.Result('new_infections_risk_group_' + str(risk_group) + '_' + sex, dtype=int, auto_plot=False),
                    ]
        self.define_results(*[res for res in results if res.name not in self.results])
        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        t = self.tally()
        self.update_sti_results(t)
        infected, new_inf, on_art, diagnosed = t.vals(self.infected), t.vals(self.new_infected()), t.vals(self.on_art), t.vals(self.diagnosed)
        res['prevalence'][ti] = t.div(t.count(infected), t.count(t.vals(ppl.alive))) # For the whole population
        res['on_art'][ti] = t.count(on_art)
        res['new_deaths'][ti] = np.count_nonzero(self.ti_dead.raw == ti)
        res['new_diagnoses'][ti] = t.count(t.vals(self.ti_diagnosed == ti))
        res['new_agents_on_art'][ti] = t.count(t.vals(self.art.ti_art == ti) & on_art) if self.art is not None else 0 # Gated on on_art, since ti_art may also hold a scheduled future start
        if self.include_mtct:
            pregnant = t.vals(ppl.pregnancy.P)
            res['n_on_art_pregnant'][ti] = t.count(on_art & pregnant)
            res['p_diagnosed_pregnant'][ti] = t.div(t.count(diagnosed & pregnant & infected), t.count(infected & pregnant))
        res['p_on_art'][ti] = t.div(res['on_art'][ti], res['infected'][ti])

        # By age, sex work, and risk group
        age = np.asarray(ppl.age.values)
        adult = (age >= 15) & (age < 50)
        res['prevalence_15_49'][ti] = t.div(t.count(infected & adult), t.count(adult))
        if 'structuredsexual' in self.sim.networks.keys():
            nw = self.sim.networks.structuredsexual
            for group, key in [(t.vals(nw.fsw), 'sw'), (t.vals(nw.client), 'client')]:
                if t.count(group):
                    res[f'prevalence_{key}'][ti] = t.div(t.count(infected & group), t.count(group))
                    res[f'new_infections_{key}'][ti] = t.count(new_inf & group)
                    res[f'new_infections_not_{key}'][ti] = t.count(new_inf & ~group)
            risk_group = np.asarray(nw.risk_group.values)
            for rg in range(nw.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    in_group = (risk_group == rg) & t.vals(ppl[sex])
                    res[f'prevalence_risk_group_{rg}_{sex}'][ti] = t.div(t.count(infected & in_group), t.count(in_group))
                    res[f'new_infections_risk_group_{rg}_{sex}'][ti] = t.count(new_inf & in_group)
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        res['cum_deaths'][:] = np.cumsum(res['new_deaths'][:])
        res['cum_diagnoses'][:] = np.cumsum(res['new_diagnoses'][:])
        return


def _art_attr(name):
    def get(self):
        art = self.art
        if art is None:
            errormsg = f'hiv.{name} is a state of treatment, held by sti.ART, and there is no ART in this sim'
            raise AttributeError(errormsg)
        return getattr(art, name)
    get.__doc__ = f'The ART module\'s `{name}` (see sti.ART)'
    return property(get)

for _name in art_attrs:
    setattr(HIV, _name, _art_attr(_name))
