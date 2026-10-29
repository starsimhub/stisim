"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
from stisim.diseases.sti import STIResults, init_prev_by_group, sexual_beta, hetero_networks, default_age_bins, default_sex_keys

__all__ = ['HIV']


class HIV(ss.Infection, STIResults):
    """
    HIV: its natural history, CD4 count, care seeking, diagnosis, treatment, PrEP, and circumcision

    People infected are in acute infection for `dur_acute`, then latent for `dur_latent`, then have a falling
    CD4 count for `dur_falling`, after which they die of AIDS (with probability `p_aids_death`); these times
    are drawn at infection. On ART, they are virally suppressed (`on_effective_art`) or not
    (`on_nonsuppressive_art`); after stopping (`art_discontinued`), they die of AIDS once their CD4 count falls
    to zero, a time worked out from their CD4 count when they stop. Throughout, they die of HIV-related illness
    at a rate by CD4 count (`cd4_death_rates`), lower on ART (`rel_art_mortality_*`). Infection is by sexual
    transmission (`beta_m2f`), and from mother to child in pregnancy (`beta_m2c`) and breastfeeding (`beta_breastfeed`).

    The CD4 count falls in acute and late infection (by the share of the stage elapsed), rises on ART, and falls
    after it. Infectiousness is by stage and ART, and susceptibility by PrEP, circumcision, and age. The ART
    intervention calls `start_art()` and `stop_art()`.

    States: `susceptible`, `infected`, `acute`, `latent`, `falling`, `on_effective_art`, `on_nonsuppressive_art`,
    and `art_discontinued`; `on_art` and `aids` are conditions over them.

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
            effective_art_efficacy = 0.99, # Reduction in transmission on effective (virally suppressive) ART: PARTNER/PARTNER2/Opposites Attract ("U=U"), less a small residual
            nonsupp_art_efficacy = 0.35, # Reduction in transmission on non-suppressive ART: Quinn (Rakai 2000), applied to an assumed non-suppressed viral load distribution
            rate_art_dropout = ss.peryear(0), # Rate of stopping ART early (e.g. from poor adherence), besides at the end of dur_on_art; none by default
            time_to_art_efficacy = ss.months(6), # Time to reach full ART efficacy: a linear increase in efficacy
            p_effective_art = ss.bernoulli(p=1.0), # Probability that a newly-initiated agent achieves viral suppression
            dur_on_art = ss.lognorm_ex(ss.years(3), ss.years(1.5)), # Base ART duration (scaled by rel_dur_on_art and trend)
            rel_dur_on_art = 1.0, # Calibratable scalar that scales ART duration
            dur_on_art_trend = None, # Optional time-varying scale, e.g. sc.objdict(years=np.array([2004, 2015, 2030]), vals=np.array([0.5, 1.0, 1.5]))
            dur_post_art = ss.normal(), # Time from stopping ART until the CD4 count reaches zero; set for each person as they stop (see stop_art()), in years
            dur_post_art_scale_factor = 0.1,
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
            ss.BoolState('on_effective_art'), # On ART and virally suppressed
            ss.BoolState('on_nonsuppressive_art'), # On ART but not suppressed
            ss.BoolState('art_discontinued'), # Off ART, having been on it
            ss.TimeArr('ti_latent'), # When each stage starts, drawn at infection
            ss.TimeArr('ti_falling'),
            ss.TimeArr('ti_zero'), # When the CD4 count reaches zero, at the end of late infection or after stopping ART: an AIDS death
            ss.FloatArr('ti_dead'), # Time of death of HIV
            ss.FloatArr('dur_post_art'), # The time from stopping ART until the CD4 count reaches zero, in steps

            # Care and treatment
            ss.FloatArr('baseline_care_seeking'),
            ss.FloatArr('care_seeking'),
            ss.BoolState('art_naive', default=True), # Never initiated ART
            ss.FloatArr('ti_art'), # Time of starting ART, or of a scheduled start (see sti.HIVTest)
            ss.TimeArr('ti_stop_art'),

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
            ss.TimeArr('ti_prep_stop'), # Scheduled course expiry (renewal/lapse), like ti_stop_art

            # Circumcision (VMMC)
            ss.BoolState('circumcised'),
            ss.BoolState('circ_traditional', default=False), # True=traditional (non-program), False=program (medical) VMMC
            ss.FloatArr('ti_circumcised'),
        )
        self.u_death = ss.random(name='u_death') # One uniform a step for each person infected, for death
        self.u_dropout = ss.random(name='u_dropout')
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

    def death_hazard(self, uids):
        """ The annual hazard of death from HIV-related illness of each agent, by CD4 count, lower on ART (not the AIDS death when the CD4 count reaches zero, at a set time), e.g. for MIGHTI's competing risks of death """
        uids = ss.uids(uids)
        h = np.zeros(len(uids))
        off, on = self.infected[uids] & ~self.on_art[uids], self.on_art[uids]
        h[off] = self.off_art_death_rate(uids[off])
        h[on] = self.on_art_death_rate(uids[on])
        return h

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
        """ CD4 counts and care seeking; and the infections at the start, which began before it (`dist_ti_init_infected`) """
        self.init_cd4()
        self.init_care_seeking()
        ss.Disease.init_post(self) # Not ss.Infection's, which seeds the infections at the start of the run
        ppl = self.sim.people
        seeds = self.pars.init_prev.filter(ppl.auids[ppl.age.raw[ppl.auids] >= 0]) # Not the fetuses of pregnancies under way
        self.set_prognoses(seeds, ti=self.pars.dist_ti_init_infected.rvs(seeds).astype(int))
        diagnosed = self.pars.init_diagnosed.filter(seeds)
        self.diagnosed[diagnosed] = True
        self.ti_diagnosed[diagnosed] = 0
        self.ti_art[diagnosed] = 0 # Scheduled to start ART at once, so an ART intervention picks them up
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
        ti = self.ti
        share = lambda start, end: np.clip(np.divide(ti - start, end - start, out=np.ones(len(start)), where=end > start), 0, 1) # Of a stage elapsed
        acute = self.acute.uids
        cd4_start, cd4_latent = self.cd4_start[acute], self.cd4_latent[acute]
        self.cd4[acute] = cd4_start - (cd4_start - cd4_latent)*share(self.ti_infected[acute], self.ti_latent[acute])
        latent = self.latent.uids
        self.cd4[latent] = self.cd4_latent[latent]
        falling = self.falling.uids
        self.cd4[falling] = np.maximum(1, self.cd4_latent[falling] - (self.cd4_latent[falling] - 1)*share(self.ti_falling[falling], self.ti_zero[falling]))
        on_art = self.on_art.uids
        self.cd4[on_art] = self.cd4_increase(on_art)
        stopped = self.art_discontinued.uids
        self.cd4[stopped] = np.maximum(1, self.cd4_postart[stopped] - (self.cd4_postart[stopped] - 1)*share(self.ti_stop_art[stopped], self.ti_zero[stopped]))
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
        art_uids = self.on_art.uids
        if len(art_uids):
            self.rel_trans[art_uids] *= self.art_factor(art_uids)
        return

    def art_factor(self, uids):
        """ The multiplier of ART on the infectiousness of some agents on it, which reaches its full efficacy by a linear rise over time_to_art_efficacy """
        p = self.pars
        full_eff = np.where(self.on_effective_art[uids], p.effective_art_efficacy, p.nonsupp_art_efficacy)
        return 1 - full_eff*np.minimum((self.ti - self.ti_art[uids])/(p.time_to_art_efficacy/self.dt), 1)

    def infectiousness(self, uids=None):
        """ How infectious each agent is, all told, relative to untreated latent infection (their rel_trans) """
        return self.rel_trans[self.sim.people.auids if uids is None else ss.uids(uids)]

    def step_state(self):
        """ CD4 counts and care seeking; the natural history; death; then each agent's susceptibility and infectiousness """
        ti = self.ti
        self.init_cd4()
        self.init_care_seeking()
        if self.include_mtct: # Pregnant women seek care more, which makes them much less likely to stop treatment
            pregnant = self.sim.demographics.pregnancy.P
            self.care_seeking[pregnant] = self.baseline_care_seeking[pregnant] * self.pars.maternal_care_scale
            self.care_seeking[~pregnant] = self.baseline_care_seeking[~pregnant]

        # The natural history, untreated: ART takes people out of it (see start_art())
        latent = (self.acute & self.ti_latent.due()).uids
        self.acute[latent] = False
        self.latent[latent] = True
        falling = (self.latent & self.ti_falling.due()).uids
        self.latent[falling] = False
        self.falling[falling] = True
        self.update_cd4()

        # Stopping ART early
        rate = self.pars.rate_art_dropout
        if rate.value > 0:
            on_art = self.on_art.uids
            self.stop_art(on_art[self.u_dropout.rvs(on_art) < rate.to_prob(self.dt)])

        # Death: of AIDS once the CD4 count reaches zero, and at any time at a rate by CD4 count, lower on ART
        infected = self.infected.uids
        on_art = self.on_art[infected]
        p_death = np.empty(len(infected))
        p_death[~on_art] = ss.peryear(self.off_art_death_rate(infected[~on_art])).to_prob(self.dt)
        p_death[on_art] = ss.peryear(self.on_art_death_rate(infected[on_art])).to_prob(self.dt)
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
        for state in [self.infected, self.acute, self.latent, self.falling, self.on_effective_art, self.on_nonsuppressive_art, self.art_discontinued, self.on_prep, self.susceptible]:
            state[uids] = False
        return

    def new_infected(self):
        """ Who was infected in this step (not those at the start, whose infections began before it) """
        return self.ti_infected == self.ti


    # Treatment, called by the interventions
    def start_art(self, uids, p_effective_art=None):
        """
        Put agents on ART: suppressed, or not; it takes them out of the natural history

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
        for state in [self.acute, self.latent, self.falling, self.art_discontinued]:
            state[uids] = False
        self.on_effective_art[uids] = False
        self.on_effective_art[effective] = True
        self.on_nonsuppressive_art[uids] = True
        self.on_nonsuppressive_art[effective] = False
        later = uids[self.ti_zero[uids] > ti] # Those due to die of AIDS now still do
        self.ti_zero[later] = np.nan
        self.ti_latent[uids] = self.ti_falling[uids] = np.nan

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
        """ Take agents off ART (by default, those whose time on it is up): they die of AIDS once their CD4 count falls to zero, a time that depends on it """
        ti = self.ti
        uids = (self.on_art & self.ti_stop_art.due()).uids if uids is None else ss.uids(uids)
        on = uids[self.on_art[uids]]
        factor = self.art_factor(on)
        self.rel_trans[on] = np.divide(self.rel_trans[on], factor, out=np.ones(len(on)), where=factor > 0) # Off ART from this step's transmission on, as they are put on it (see update_transmission())
        self.on_effective_art[uids] = False
        self.on_nonsuppressive_art[uids] = False
        self.art_discontinued[uids] = True
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
        return

    def set_suppression(self, uids, effective):
        """ Make agents on ART virally suppressed (effective=True), or not """
        uids = ss.uids(uids)
        uids = uids[self.on_art[uids]]
        self.on_effective_art[uids] = effective
        self.on_nonsuppressive_art[uids] = not effective
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
        res['new_agents_on_art'][ti] = t.count(t.vals(self.ti_art == ti) & on_art) # Gated on on_art, since ti_art may also hold a scheduled future start
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
