"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim.utils as ut
from stisim.diseases.sti import STIResults, init_prev_by_group, sexual_beta, hetero_networks, default_age_bins, default_sex_keys

__all__ = ['HIV']


class HIV(ss.Disease, STIResults):
    """
    HIV, with its natural history and care states declared, and its CD4 count, care seeking, diagnosis, PrEP, and circumcision as code

    Declared: people infected are in acute infection for `dur_acute`, then latent for `dur_latent`, then have a
    falling CD4 count for `dur_falling`, after which they die of AIDS (with probability `p_aids_death`). On ART,
    they are virally suppressed (`on_effective_art`) or not (`on_nonsuppressive_art`); after stopping
    (`art_discontinued`), they die of AIDS once their CD4 count falls to zero, a time worked out from their CD4
    count when they stop. Throughout, they die of HIV-related illness at a rate by CD4 count (`cd4_death_rates`),
    lower on ART (`rel_art_mortality_*`). Infection is by sexual transmission (`beta_m2f`), and from mother to
    child in pregnancy (`beta_m2c`) and breastfeeding (`beta_breastfeed`).

    As code: the CD4 count, which declines in acute and falling infection (by the share of the stage elapsed),
    rises on ART, and falls after it; infectiousness by stage and ART, and susceptibility by PrEP, circumcision,
    and age; care seeking; diagnosis; and starting and stopping ART (`start_art()`, `stop_art()`), which the
    ART intervention calls.

    States: `susceptible`, `acute`, `latent`, `falling`, `on_effective_art`, `on_nonsuppressive_art`, and
    `art_discontinued`; `infected`, `on_art`, and `aids` are conditions over them.

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
                 sex_keys=default_sex_keys, include_care=True, **kwargs):
        values = sc.mergedicts(dict(init_prev=0.05), pars, kwargs)
        code_pars = self.code_defaults()
        given = {key:values.pop(key) for key in code_pars if key in values}
        self.init_prev_data = init_prev_data
        if init_prev_data is not None:
            values['init_prev'] = self.init_prev_by_group
        super().__init__(name=name, **sc.mergedicts(dict(dt='month'), values)) # Monthly unless dt is given

        # Transmission: a probability per act on the sexual networks, and rates per month from mother to child
        beta_m2f = ss.par(0.05, doc='Probability of transmission per act from a man to a woman')
        rel_beta_f2m = ss.par(0.5, doc='Transmission from a woman to a man, relative to from a man to a woman')
        beta_m2m = ss.par(0.0, doc='Probability of transmission per act between men')
        beta_m2c = ss.par(0.025, doc='Rate of transmission from mother to child in pregnancy, per month (about 20% over a pregnancy)')
        beta_breastfeed = ss.par(0.005, doc='Rate of transmission from mother to child by breastfeeding, per month (about 14% over 12 months without ART)')

        # Natural history
        dur_acute = ss.par(ss.lognorm_ex(ss.months(1.7), ss.months(1)), doc='Duration of acute infection') # Bellan 2015: central estimate 1.7 months
        dur_latent = ss.par(ss.lognorm_ex(ss.years(10), ss.years(3)), doc='Duration of latent, untreated infection')
        dur_falling = ss.par(ss.lognorm_ex(ss.years(3), ss.years(1)), doc='Duration of late-stage infection, as the CD4 count falls to zero')
        p_aids_death = ss.par(1.0, doc='Probability of dying of AIDS when the CD4 count reaches zero')

        # Treatment: how much ART reduces transmission, and stopping it early (e.g. from poor adherence), besides at the end of dur_on_art
        effective_art_efficacy = ss.par(0.99, doc='Reduction in transmission on effective (virally suppressive) ART') # PARTNER/PARTNER2/Opposites Attract ("U=U"), less a small residual
        nonsupp_art_efficacy = ss.par(0.35, doc='Reduction in transmission on non-suppressive ART') # Quinn (Rakai 2000), applied to an assumed non-suppressed viral load distribution
        rate_art_dropout = ss.par(ss.peryear(0), doc='Rate of stopping ART early; none by default')

        # States
        self.define_states(
            ss.State('susceptible'),
            ss.State('acute', keep_time=True), # The times are kept for the CD4 count (see step_state())
            ss.State('latent'),
            ss.State('falling', keep_time=True),
            ss.State('on_effective_art'), # On ART and virally suppressed
            ss.State('on_nonsuppressive_art'), # On ART but not suppressed
            ss.State('art_discontinued', keep_time=True), # Off ART, having been on it
        )
        acute, latent, falling = self.acute, self.latent, self.falling
        on_art = self.on_effective_art + self.on_nonsuppressive_art
        off_art = acute + latent + falling + self.art_discontinued
        infectious = off_art + (1 - effective_art_efficacy)*self.on_effective_art + (1 - nonsupp_art_efficacy)*self.on_nonsuppressive_art # ART at full efficacy; see update_transmission() for until then

        # Transmission
        sexual = sexual_beta(beta_m2f, beta_m2f*rel_beta_f2m, beta_m2m, hetero=self.hetero_networks)
        self.define_transitions(
            infection = self.susceptible >> acute @ sexual*infectious,
            mtct_prenatal = self.susceptible >> acute @ beta_m2c*dict(maternalnet=ss.permonth(1), default=0)*infectious,
            mtct_postnatal = self.susceptible >> acute @ beta_breastfeed*dict(breastfeedingnet=ss.permonth(1), default=0)*infectious,
        )

        # Natural history, and death: of AIDS once the CD4 count reaches zero, and at any time at a rate by CD4 count (see death_rate())
        self.define_transitions(
            latency = acute >> latent @ dur_acute,
            progression = latent >> falling @ dur_latent,
            aids_death = ss.Transition(falling, [ss.null, ss.stay], p=[p_aids_death, None], dur=dur_falling),
            post_art_death = ss.Transition(self.art_discontinued, [ss.null, ss.stay], p=[p_aids_death, None], dur=ss.constant(self.post_art_dwell)), # The time drawn as they stop ART (see post_art_dwell())
            hiv_death = ss.Transition(off_art, ss.null, rate=self.death_rate),
            art_death = ss.Transition(on_art, ss.null, rate=self.art_death_rate),
        )

        # Stopping ART early (e.g. ss.modify(hiv.pars.rate_art_dropout, set=...) for poor adherence); stop_art() takes people off it otherwise
        self.define_transitions(
            art_dropout = ss.Transition(on_art, self.art_discontinued, rate=rate_art_dropout),
        )

        # Results of the flows
        self.define_results(
            new_deaths = self.aids_death + self.post_art_death + self.hiv_death + self.art_death,
            new_infections_sex = 1*self.infection,
            new_infections_prenatal = 1*self.mtct_prenatal,
            new_infections_postnatal = 1*self.mtct_postnatal,
            new_infections_mtct = self.mtct_prenatal + self.mtct_postnatal,
        )

        # As code: the CD4 count, care, diagnosis, PrEP, and circumcision
        self.define_states(
            ss.FloatArr('ti_infected'), # Time of infection, as it was
            ss.FloatArr('ti_dead'), # Time of death of HIV
            ss.FloatArr('ti_zero'), # Time the CD4 count reaches zero after stopping ART, which is generally an AIDS death
            ss.FloatArr('dur_post_art'), # The time from stopping ART until then, in steps

            # Care and treatment
            ss.FloatArr('baseline_care_seeking'),
            ss.FloatArr('care_seeking'),
            ss.BoolState('art_naive', default=True), # Never initiated ART
            ss.FloatArr('ti_art'), # Time of starting ART, or of a scheduled start (see sti.HIVTest)
            ss.FloatArr('ti_stop_art'),

            # CD4 count
            ss.FloatArr('cd4'), # Current CD4 count
            ss.FloatArr('cd4_start'), # Initial CD4 count for each agent before an infection
            ss.FloatArr('cd4_preart'), # CD4 immediately before initiating ART
            ss.FloatArr('cd4_latent'), # CD4 count during latent infection
            ss.FloatArr('cd4_nadir'), # Lowest CD4
            ss.FloatArr('cd4_potential'), # Potential CD4 count if continually treated
            ss.FloatArr('cd4_postart'), # CD4 after stopping ART

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
            ss.FloatArr('ti_prep_stop'), # Scheduled course expiry (renewal/lapse), like ti_stop_art

            # Circumcision (VMMC)
            ss.BoolState('circumcised'),
            ss.BoolState('circ_traditional', default=False), # True=traditional (non-program), False=program (medical) VMMC
            ss.FloatArr('ti_circumcised'),
        )
        self.define_pars(**code_pars)
        self.pars.update(given) # Including those with bounds, which are kept for fitting
        self.set_strata([15, 50], age_bins, sex_keys)
        self.include_care = include_care
        return

    @staticmethod
    def code_defaults():
        """ The parameters that the disease's own code reads, rather than its declared transitions """
        return dict(
            # Death from HIV-related illness, by CD4 count, off ART
            cd4_death_bins = np.array([1000, 500, 350, 200, 50, 0]), # Off-ART CD4-stratified mortality (descending)
            cd4_death_rates = np.array([0.003, 0.003, 0.005, 0.01, 0.05, 0.30]), # Annual rate per bin
            rel_death = 1.0, # Scales all HIV death rates (off- and on-ART)
            rel_death_f = 0.74, # Additional mortality multiplier for females, on- and off-ART alike

            # Mortality on ART, relative to the off-ART CD4-based rate, so that being on ART is never more dangerous than being
            # off it at the same CD4 count (as long as the largest of these times the largest age and duration multipliers is
            # at most 1, which init_pre() checks)
            rel_art_mortality_effective = 0.25, # Fraction of off-ART mortality retained on effective (suppressive) ART
            rel_art_mortality_unsupp_m = 0.7, # ... on non-suppressive ART, males
            rel_art_mortality_unsupp_f = 0.35, # ... females: the non-suppressive/effective ratio is 2x higher for men than for women
            art_death_age = [(0, 25, 1.0), (25, 35, 1.10), (35, 45, 1.21), (45, 125, 1.32)], # (age_lo, age_hi, mult)
            art_death_dur = None, # (dur_lo, dur_hi, mult) list, days since ti_art; None = no-op

            # CD4 count
            cd4_start = ss.normal(loc=800, scale=50),
            cd4_latent = ss.normal(loc=500, scale=50),
            art_cd4_growth = 0.1, # Unitless parameter of how CD4 reconstitutes after starting ART, in a logistic growth function
            art_cd4_pars = dict(cd4_max=1000, cd4_healthy=500),

            # Infectiousness and susceptibility
            eff_condom = 0.9, # Read by the sexual networks (see sti.BaseNetwork.net_beta())
            rel_trans_acute = ss.normal(loc=5.3, scale=0.5), # Bellan 2015: RH=5.3, EHM≈7.3
            rel_trans_falling = ss.normal(loc=8, scale=0.5), # Increase in late infection (a CD4 count below 200)
            rel_sus_age = None, # Age/sex-stratified rel_sus multipliers: (age_lo, age_hi, sex, multiplier) tuples; sex is 'f', 'm', or None for both. E.g. [(15, 25, 'f', 1.7), (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)]
            eff_circ = 0.6, # Program (medical) VMMC: reduction in HIV acquisition susceptibility
            eff_circ_traditional = 0.6, # Traditional (non-program) circumcision: reduction in HIV acquisition susceptibility

            # Initial conditions
            rel_init_prev = 1, # Scales init_prev_data
            init_diagnosed = ss.bernoulli(p=0),
            dist_ti_init_infected = ss.uniform(low=-10*12, high=-5), # When the infections at the start began, in steps

            # Care seeking
            care_seeking = ss.normal(loc=1, scale=0.5), # Relative care-seeking behavior (wide spread for low-care-seeking extremes)
            maternal_care_scale = 2, # Factor for scaling up care-seeking behavior during pregnancy

            # Treatment
            time_to_art_efficacy = ss.months(6), # Time to reach full ART efficacy: a linear increase in efficacy
            p_effective_art = ss.bernoulli(p=1.0), # Probability that a newly-initiated agent achieves viral suppression
            dur_on_art = ss.lognorm_ex(ss.years(3), ss.years(1.5)), # Base ART duration (scaled by rel_dur_on_art and trend)
            rel_dur_on_art = 1.0, # Calibratable scalar that scales ART duration
            dur_on_art_trend = None, # Optional time-varying scale, e.g. sc.objdict(years=np.array([2004, 2015, 2030]), vals=np.array([0.5, 1.0, 1.5]))
            dur_post_art = ss.normal(), # Time from stopping ART until the CD4 count reaches zero; set for each person as they stop (see post_art_dwell()), in years
            dur_post_art_scale_factor = 0.1,
        )

    # Conditions over the states
    @property
    def infected(self):
        """ Infected: in any state but susceptible """
        return self.acute | self.latent | self.falling | self.on_art | self.art_discontinued

    @property
    def infectious(self):
        return self.infected

    @property
    def symptomatic(self):
        return self.infectious

    @property
    def treatable(self):
        return self.infected

    @property
    def on_art(self):
        """ On ART, suppressed or not """
        return self.on_effective_art | self.on_nonsuppressive_art

    @property
    def aids(self):
        return self.cd4 < 200

    @property
    def include_mtct(self):
        return 'pregnancy' in self.sim.demographics

    def init_prev_by_group(self, sim, uids):
        """ The initial prevalence of each agent, from init_prev_data """
        return init_prev_by_group(sim, uids, self.init_prev_data, self.pars.rel_init_prev)

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

    # Death, as the rates of the declared transitions
    def off_art_death_rate(self, uids):
        """ The annual rate of death from HIV-related illness of each agent off ART, by CD4 count, sex, and rel_death """
        rate = self.pars.cd4_death_rates[np.digitize(self.cd4[uids], self.pars.cd4_death_bins)] * self.pars.rel_death
        rate[~self.sim.people.male[uids]] *= self.pars.rel_death_f
        return rate

    def on_art_death_rate(self, uids):
        """
        The annual rate of death of each agent on ART: the off-ART rate at the same CD4 count, times rel_art_mortality
        (by suppression, and by sex if not suppressed) and the multipliers by age (art_death_age) and time on ART (art_death_dur)
        """
        male = self.sim.people.male[uids]
        effective = self.on_effective_art[uids]
        age = self.sim.people.age[uids]
        rate = self.off_art_death_rate(uids)
        rate[effective] *= self.pars.rel_art_mortality_effective
        rate[~effective & male] *= self.pars.rel_art_mortality_unsupp_m
        rate[~effective & ~male] *= self.pars.rel_art_mortality_unsupp_f

        age_mult = np.ones(len(uids))
        for age_lo, age_hi, mult in self.pars.art_death_age:
            age_mult[(age >= age_lo) & (age < age_hi)] *= mult
        rate = rate * age_mult

        if self.pars.art_death_dur is not None:
            dur_days = (self.ti - self.ti_art[uids]) * self.dt.days
            dur_mult = np.ones(len(uids))
            for dur_lo, dur_hi, mult in self.pars.art_death_dur:
                dur_mult[(dur_days >= dur_lo) & (dur_days < dur_hi)] *= mult
            rate = rate * dur_mult
        return rate

    def death_rate(self, sim, uids):
        """ The rate of the declared transition hiv_death """
        return ss.peryear(self.off_art_death_rate(uids))

    def art_death_rate(self, sim, uids):
        """ The rate of the declared transition art_death """
        return ss.peryear(self.on_art_death_rate(uids))

    def post_art_dwell(self, sim, uids):
        """
        As people stop ART (by stop_art() or art_dropout): the time until their CD4 count falls to zero, which depends on it, and when they die of AIDS (post_art_death)
        """
        ti = self.ti
        self.ti_stop_art[uids] = ti
        self.cd4_postart[uids] = self.cd4[uids]
        dur_mean = np.log(self.cd4_preart[uids])*self.cd4[uids]/self.cd4_potential[uids]
        dur_scale = dur_mean*self.pars.dur_post_art_scale_factor
        self.pars.dur_post_art.set(loc=ss.years(dur_mean), scale=ss.years(np.maximum(dur_scale, 1e-3))) # Not zero
        dur_post_art = self.pars.dur_post_art.rvs(uids)
        if np.isnan(dur_post_art).any():
            bad = uids[np.isnan(dur_post_art)]
            errormsg = f'Invalid post-ART duration for {bad}: CD4 {self.cd4[bad]}, before ART {self.cd4_preart[bad]}, potential {self.cd4_potential[bad]}'
            raise ValueError(errormsg)
        self.dur_post_art[uids] = np.maximum(dur_post_art.astype(int), 0)
        self.ti_zero[uids] = ti + self.dur_post_art[uids]
        return self.dt*self.dur_post_art[uids]

    def init_pre(self, sim):
        super().init_pre(sim)
        # On-ART mortality is anchored to the off-ART CD4-based rate so it can never exceed it by construction, as long as
        # the largest combination of multipliers stays <= 1. Catch a violation here rather than let a calibration silently
        # invert the treatment effect.
        p = self.pars
        age_mult = max(mult for *_, mult in p.art_death_age) if p.art_death_age else 1.0
        dur_mult = max(mult for *_, mult in p.art_death_dur) if p.art_death_dur else 1.0
        worst = max(p.rel_art_mortality_effective, p.rel_art_mortality_unsupp_m, p.rel_art_mortality_unsupp_f) * age_mult * dur_mult
        if worst > 1:
            errormsg = (f'On-ART mortality can exceed off-ART mortality at the same CD4 count: max(rel_art_mortality_effective, '
                        f'rel_art_mortality_unsupp_m, rel_art_mortality_unsupp_f) * max(art_death_age mult) * max(art_death_dur mult) '
                        f'= {worst:.3f} > 1. Reduce these pars so the product stays <= 1.')
            raise ValueError(errormsg)
        return

    def init_post(self):
        """ CD4 counts and care seeking; and the infections at the start, which began before it """
        self.init_cd4()
        self.init_care_seeking()
        super().init_post()
        seeds = self.acute.uids
        self.backdate_seeds(seeds)
        self.set_infection(seeds)
        diagnosed = self.pars.init_diagnosed.filter(seeds)
        self.diagnosed[diagnosed] = True
        self.ti_diagnosed[diagnosed] = 0
        self.ti_art[diagnosed] = 0 # Scheduled to start ART at once, so an ART intervention picks them up
        return

    def backdate_seeds(self, uids):
        """
        The infections at the start began dist_ti_init_infected before it: each is put in the stage it has reached,
        partway through it, with the stage's durations drawn as for any infection
        """
        if not len(uids):
            return
        eng = self.engine
        trs = self.graph.transitions
        now = eng.ti
        began = now + self.pars.dist_ti_init_infected.rvs(uids).astype(int) # In steps, before the start
        steps = lambda tr, state: eng.whole_steps(state, uids, eng.sample_dwell(tr, uids))
        to_latent = began + steps(trs.latency, 'acute')
        to_falling = to_latent + steps(trs.progression, 'latent')
        to_zero = to_falling + steps(trs.aids_death, 'falling')
        self.ti_infected[uids] = began
        for state, start, end in [('acute', began, to_latent), ('latent', to_latent, to_falling), ('falling', to_falling, to_zero)]:
            here = (start <= now) & ((end > now) | (state == 'falling')) # Those whose course ended before the start are due to die in its first step
            if here.any():
                self.move(uids[here], self[state])
                eng.backdate(state, uids[here], entered=start[here], exit=end[here]) # [PROPOSAL ST-E] ss.move(..., since=)
        return

    # CD4 count and care seeking
    def init_cd4(self):
        """ Set the CD4 counts of new agents """
        uids = ss.uids(np.isnan(self.cd4_start.raw).nonzero()[-1])
        self.cd4_start[uids] = self.pars.cd4_start.rvs(uids)
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        return

    def init_care_seeking(self):
        """ Set the care-seeking behavior of new agents """
        uids = ss.uids(self.care_seeking.isnan)
        self.care_seeking[uids] = self.pars.care_seeking.rvs(uids)
        self.care_seeking[uids] = np.maximum(self.care_seeking[uids], 0.01) # Floor to prevent division by zero
        self.baseline_care_seeking[uids] = sc.dcp(self.care_seeking[uids]) # Copy it so pregnancy can modify it
        return

    def progress(self, state, uids):
        """ The share of their time in a declared state that agents have spent, from 0 to 1 """
        return self.engine.clock_values(self.graph.states[state].progress, uids)

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
        dur_art = self.ti - self.ti_art[uids]
        cd4_preart = self.cd4_preart[uids]
        cd4_total_gain = self.cd4_potential[uids] - cd4_preart
        return 2*cd4_total_gain/(1 + np.exp(-dur_art*self.pars.art_cd4_growth)) - cd4_total_gain + cd4_preart # Concave logistic

    def update_cd4(self):
        """ CD4 counts: falling from the start of infection to the latent level in acute infection, and from it to 1 in late infection; rising on ART; and falling to 1 after it """
        acute = self.acute.uids
        cd4_start, cd4_latent = self.cd4_start[acute], self.cd4_latent[acute]
        self.cd4[acute] = cd4_start - (cd4_start - cd4_latent)*self.progress('acute', acute)
        latent = self.latent.uids
        self.cd4[latent] = self.cd4_latent[latent]
        falling = self.falling.uids
        self.cd4[falling] = np.maximum(1, self.cd4_latent[falling] - (self.cd4_latent[falling] - 1)*self.progress('falling', falling))
        on_art = self.on_art.uids
        self.cd4[on_art] = self.cd4_increase(on_art)
        stopped = self.art_discontinued.uids
        self.cd4[stopped] = np.maximum(1, self.cd4_postart[stopped] - (self.cd4_postart[stopped] - 1)*self.progress('art_discontinued', stopped))
        untreated = (self.infected & ~self.on_art).uids
        self.cd4_nadir[untreated] = np.minimum(self.cd4_nadir[untreated], self.cd4[untreated])
        return

    def update_transmission(self):
        """
        Each agent's susceptibility (age, PrEP, circumcision) and infectiousness (acute and late infection; for ART, see art_ramp()), reset each step; connectors change them after this
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

        # Infectiousness in acute and late infection
        self.rel_trans[self.acute] *= p.rel_trans_acute.rvs(self.acute.uids)
        self.rel_trans[self.aids] *= p.rel_trans_falling.rvs(self.aids.uids)

        return

    def art_ramp(self):
        """
        On ART: the infection transitions count people at ART's full efficacy, which they reach by a linear rise over time_to_art_efficacy; until
        then, they are more infectious than that, by this ratio (with an efficacy of 1, there is nothing to scale, so they do not transmit from
        the start). Applied just before transmission, after the interventions that start ART in the step
        """
        art_uids = self.on_art.uids
        if len(art_uids):
            full_eff = np.where(self.on_effective_art[art_uids], self.value('effective_art_efficacy'), self.value('nonsupp_art_efficacy'))
            steps_on_art = self.ti - self.ti_art[art_uids]
            steps_to_full = self.pars.time_to_art_efficacy / self.dt
            efficacy = full_eff*np.minimum(steps_on_art/steps_to_full, 1)
            self.rel_trans[art_uids] *= np.divide(1 - efficacy, 1 - full_eff, out=np.ones(len(art_uids)), where=full_eff < 1)
        return

    def infectiousness(self, uids=None):
        """
        How infectious each agent is, all told, relative to untreated infection: their rel_trans, times the reduction of ART at its
        full efficacy, which the infection transitions apply (so rel_trans on ART is relative to that; see update_transmission())
        """
        uids = self.sim.people.auids if uids is None else ss.uids(uids)
        full = np.where(self.on_effective_art[uids], 1 - self.value('effective_art_efficacy'), 1.0)
        full = np.where(self.on_nonsuppressive_art[uids], 1 - self.value('nonsupp_art_efficacy'), full)
        return self.rel_trans[uids]*full

    def step_state(self):
        """ CD4 counts and care seeking; the declared transitions; then each agent's susceptibility and infectiousness """
        self.init_cd4()
        self.init_care_seeking()
        if self.include_mtct: # Pregnant women seek care more, which makes them much less likely to stop treatment
            pregnant = self.sim.demographics.pregnancy.P
            self.care_seeking[pregnant] = self.baseline_care_seeking[pregnant] * self.pars.maternal_care_scale
            self.care_seeking[~pregnant] = self.baseline_care_seeking[~pregnant]
        self.update_cd4()
        super().step_state()
        died = self.aids_death.uids | self.post_art_death.uids | self.hiv_death.uids | self.art_death.uids
        self.ti_dead[died] = self.ti
        self.update_transmission()
        return

    def step(self):
        """ Transmission, with ART's rise to full efficacy; then the CD4 counts of those infected (a child's from its mother) """
        self.art_ramp()
        super().step()
        children = self.mtct_prenatal.uids | self.mtct_postnatal.uids
        mothers, _ = self.engine.sources(children)
        known = ~np.isnan(mothers)
        self.cd4_start[children[known]] = self.cd4_start[ss.uids(mothers[known])]
        self.cd4_nadir[children] = self.cd4_start[children]
        infections = [self[name].uids for name,tr in self.graph.transitions.items() if 'acute' in [state.name for state in tr.dst]] # However they were infected (e.g. by a transition a subclass adds)
        self.set_infection(ss.uids.concatenate(*infections))
        return

    def set_infection(self, uids):
        """ At infection: the CD4 count, and the latent level it falls to """
        self.ti_infected[uids[np.isnan(self.ti_infected[uids])]] = self.ti # Those at the start have their own
        self.cd4[uids] = self.cd4_start[uids]
        self.cd4_latent[uids] = self.pars.cd4_latent.rvs(uids)
        return

    def new_infected(self):
        """ Who was infected in this step (not those at the start, whose infections began before it) """
        return self.ti_infected == self.ti

    # Treatment, called by the interventions
    def start_art(self, uids, p_effective_art=None):
        """
        Put agents on ART: suppressed, or not

        Args:
            uids: agents starting ART treatment
            p_effective_art (float/ss.Dist, optional): probability that each agent achieves viral suppression (effective
                ART) rather than non-suppressive ART. Pass a float to override `self.pars.p_effective_art`'s probability,
                or an already-initialized `ss.Dist` (e.g. an intervention's own par) to use it directly. Defaults to
                `self.pars.p_effective_art` (default: always effective).
        """
        ti = self.ti
        p = self.pars
        newly_treated = uids[self.art_naive[uids]]
        self.art_naive[newly_treated] = False
        self.ti_art[uids] = ti
        self.cd4_preart[uids] = self.cd4[uids]

        # Whether each achieves viral suppression (effective ART) or not, e.g. due to poor adherence or drug resistance
        if isinstance(p_effective_art, ss.Dist):
            effective = p_effective_art.filter(uids)
        else:
            if p_effective_art is not None:
                p.p_effective_art.set(p=p_effective_art)
            effective = p.p_effective_art.filter(uids)
        self.move(effective, self.on_effective_art)
        self.move(uids.remove(effective), self.on_nonsuppressive_art)

        # When they will stop
        dur_on_art = p.dur_on_art.rvs(uids) * p.rel_dur_on_art
        trend = p.dur_on_art_trend
        if trend is not None:
            dur_on_art = dur_on_art * np.interp(self.t.now('year'), trend.years, trend.vals)
        self.ti_stop_art[uids] = ti + dur_on_art.astype(int)

        # The CD4 count they can reach on treatment, set at their first start: growth follows the concave part of a logistic function, and the total gain depends on the CD4 count at initiation
        if len(newly_treated):
            cd4_max = p.art_cd4_pars['cd4_max']
            cd4_healthy = p.art_cd4_pars['cd4_healthy']
            cd4_preart = self.cd4_preart[newly_treated]
            if (cd4_preart == 0).any():
                raise ValueError('CD4 count is zero at initiation of ART')
            cd4_scale_factor = (cd4_max - cd4_preart)/cd4_healthy*np.log(cd4_max/cd4_preart)
            self.cd4_potential[newly_treated] = cd4_preart + cd4_preart*cd4_scale_factor
        return

    def stop_art(self, uids=None):
        """ Take agents off ART (by default, those whose time on it is up); see post_art_dwell() for what follows """
        if uids is None:
            uids = (self.on_art & (self.ti_stop_art <= self.ti)).uids
        self.move(ss.uids(uids), self.art_discontinued)
        return

    def set_suppression(self, uids, effective):
        """ Make agents on ART virally suppressed (effective=True), or not """
        uids = ss.uids(uids)
        if effective:
            self.move(uids[self.on_nonsuppressive_art[uids]], self.on_effective_art)
        else:
            self.move(uids[self.on_effective_art[uids]], self.on_nonsuppressive_art)
        return

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
            uids = (self.on_prep & (self.ti_prep_stop <= self.ti)).uids
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
            ss.Result('infected', dtype=int, label='Number infected'),
            ss.Result('on_art', dtype=int, label='On ART'),
            ss.Result('diagnosed', dtype=int, label='Diagnosed', auto_plot=False),
            ss.Result('cum_infections', dtype=int, label='Cumulative infections', auto_plot=False),
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
        self.define_results(*[res for res in results if res.name not in self.results and res.name not in self.graph.exprs and res.name not in self.graph.outcomes])
        return

    def update_results(self):
        super().update_results()
        self.update_sti_results()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        new_inf = self.new_infected()

        res['prevalence'][ti] = ut.cond_prob(self.infected, ppl.alive) # For the whole population
        res['infected'][ti] = ut.count(self.infected)
        res['on_art'][ti] = ut.count(self.on_art)
        res['diagnosed'][ti] = ut.count(self.diagnosed)
        res['new_diagnoses'][ti] = ut.count(self.ti_diagnosed == ti)
        res['new_agents_on_art'][ti] = ut.count((self.ti_art == ti) & self.on_art) # Gated on on_art, since ti_art may also hold a scheduled future start
        if self.include_mtct:
            pregnant = ppl.pregnancy.P
            res['n_on_art_pregnant'][ti] = ut.count(self.on_art & pregnant)
            res['p_diagnosed_pregnant'][ti] = sc.safedivide(ut.count(self.diagnosed & pregnant & self.infected), ut.count(self.infected & pregnant))
        res['p_on_art'][ti] = sc.safedivide(res['on_art'][ti], res['infected'][ti])

        # By age, sex work, and risk group
        res['prevalence_15_49'][ti] = ut.cond_prob(self.infected, (ppl.age >= 15) & (ppl.age < 50))
        if 'structuredsexual' in self.sim.networks.keys():
            nw = self.sim.networks.structuredsexual
            for group, key in [(nw.fsw, 'sw'), (nw.client, 'client')]:
                if ut.count(group):
                    res[f'prevalence_{key}'][ti] = ut.cond_prob(self.infected, group)
                    res[f'new_infections_{key}'][ti] = ut.count(new_inf & group)
                    res[f'new_infections_not_{key}'][ti] = ut.count(new_inf & ~group)
            for risk_group in range(nw.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    in_group = (nw.risk_group == risk_group) & ppl[sex]
                    res[f'prevalence_risk_group_{risk_group}_{sex}'][ti] = ut.cond_prob(self.infected, in_group)
                    res[f'new_infections_risk_group_{risk_group}_{sex}'][ti] = ut.count(new_inf & in_group)
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        res['cum_infections'][:] = np.cumsum(res['new_infections'][:])
        res['cum_deaths'][:] = np.cumsum(res['new_deaths'][:])
        res['cum_diagnoses'][:] = np.cumsum(res['new_diagnoses'][:])
        return

