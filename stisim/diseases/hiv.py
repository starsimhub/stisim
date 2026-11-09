"""
Define HIV disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim.utils as ut
from stisim.diseases.sti import STIResults, init_prev_by_group, sexual_beta, hetero_networks, default_age_bins, default_sex_keys

__all__ = ['HIV']

# The parameters of sti.ART, which HIV takes and passes on to it (see sti.ART.init_pre()), so that they can be given to either
art_pars = ['effective_art_efficacy', 'nonsupp_art_efficacy', 'rate_art_dropout', 'rel_art_mortality_effective', 'rel_art_mortality_unsupp_m',
            'rel_art_mortality_unsupp_f', 'art_death_age', 'art_death_dur', 'art_cd4_growth', 'art_cd4_pars', 'care_seeking', 'maternal_care_scale',
            'time_to_art_efficacy', 'p_effective_art', 'dur_on_art', 'rel_dur_on_art', 'dur_on_art_trend', 'dur_post_art', 'dur_post_art_scale_factor']

# The states and methods of sti.ART that can be read through HIV (e.g. hiv.ti_art)
art_attrs = ['art_naive', 'ti_art', 'ti_stop_art', 'ti_zero', 'cd4_preart', 'cd4_potential', 'cd4_postart', 'dur_post_art', 'care_seeking',
             'baseline_care_seeking', 'start_art', 'stop_art', 'set_suppression', 'on_art_death_rate']


class HIV(ss.Disease, STIResults):
    """
    HIV, with its natural history declared, and its CD4 count, diagnosis, PrEP, and circumcision as code

    Declared: people infected are in acute infection for `dur_acute`, then latent for `dur_latent`, then have a
    falling CD4 count for `dur_falling`, after which they die of AIDS (with probability `p_aids_death`).
    Throughout, they die of HIV-related illness at a rate by CD4 count (`cd4_death_rates`). Infection is by sexual
    transmission (`beta_m2f`), and from mother to child in pregnancy (`beta_m2c`) and breastfeeding (`beta_breastfeed`).

    As code: the CD4 count, which declines in acute and falling infection (by the share of the stage elapsed);
    infectiousness by stage, and susceptibility by PrEP, circumcision, and age; and diagnosis.

    Treatment is a module of its own (`sti.ART`), which registers itself with `add_treatment()`. Treatment takes
    people out of the natural history, into the states of treatment: virally suppressed (`on_effective_art`) or
    not (`on_nonsuppressive_art`), and after stopping (`art_discontinued`), when they die of AIDS once their CD4
    count falls to zero. These states are HIV's, since a person is in one state of HIV at a time, but everything
    about them comes from the treatment: who enters and leaves them (`sti.ART.start_art()`, `stop_art()`), the
    rates of death and of stopping early in them, the time from stopping to a CD4 count of zero, the CD4 counts
    of those in them, and the multiplier of their infectiousness. HIV has no parameters of treatment of its own; those
    given to it are passed on to `sti.ART`, and treatment's states and methods can be read through it
    (`hiv.ti_art`, `hiv.start_art()`, ...). Testing refers people to treatment with `refer()`.

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
    suggested_solver = ss.solvers.abm(step='wholestep') # As agents, rates over whole steps, as STIsim's diseases took them when written as code (faster; a sim's solver=ss.solvers.abm(step='midstep') overrides it)
    hetero_networks = hetero_networks # The sexual networks between men and women, by name, which a subclass can add to

    def __init__(self, pars=None, name='hiv', init_prev_data=None, age_bins=default_age_bins,
                 sex_keys=default_sex_keys, include_care=True, **kwargs):
        values = sc.mergedicts(dict(init_prev=0.05), pars, kwargs)
        code_pars = self.code_defaults()
        self.art_given = {key:values.pop(key) for key in art_pars if key in values} # Taken by sti.ART (see its init_pre())
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
        dur_acute = ss.par(ss.lognormal(ss.months(1.7), ss.months(1)), doc='Duration of acute infection') # Bellan 2015: central estimate 1.7 months
        dur_latent = ss.par(ss.lognormal(ss.years(10), ss.years(3)), doc='Duration of latent, untreated infection')
        dur_falling = ss.par(ss.lognormal(ss.years(3), ss.years(1)), doc='Duration of late-stage infection, as the CD4 count falls to zero')
        p_aids_death = ss.par(1.0, doc='Probability of dying of AIDS when the CD4 count reaches zero')

        # States
        self.define_states(
            ss.State('susceptible'),
            ss.State('acute', keep_time=True), # The times are kept for the CD4 count (see step_state())
            ss.State('latent'),
            ss.State('falling', keep_time=True),
            ss.State('on_effective_art'), # On treatment and virally suppressed (moved into and out of by the treatment, e.g. sti.ART)
            ss.State('on_nonsuppressive_art'), # On treatment but not suppressed
            ss.State('art_discontinued', keep_time=True), # Off treatment, having been on it
        )
        acute, latent, falling = self.acute, self.latent, self.falling
        on_art = self.on_effective_art + self.on_nonsuppressive_art
        off_art = acute + latent + falling + self.art_discontinued
        infectious = off_art + on_art # Treatment's reduction of infectiousness is in rel_trans (see treat_infectiousness())

        # Transmission
        sexual = sexual_beta(beta_m2f, beta_m2f*rel_beta_f2m, beta_m2m, hetero=self.hetero_networks)
        self.define_transitions(
            infection = self.susceptible >> acute @ sexual*infectious,
            mtct_prenatal = self.susceptible >> acute @ beta_m2c*dict(prenatalnet=ss.permonth(1), default=0)*infectious,
            mtct_postnatal = self.susceptible >> acute @ beta_breastfeed*dict(breastfeedingnet=ss.permonth(1), default=0)*infectious,
        )

        # Natural history, and death: of AIDS once the CD4 count reaches zero, and at any time at a rate by CD4 count (see death_rate()), on treatment at the treatment's rate
        self.define_transitions(
            latency = acute >> latent @ dur_acute,
            progression = latent >> falling @ dur_latent,
            aids_death = ss.Transition(falling, {ss.dead: p_aids_death, ss.stay: None}, dur=dur_falling),
            post_art_death = ss.Transition(self.art_discontinued, {ss.dead: p_aids_death, ss.stay: None}, dur=ss.constant(self.post_art_dwell)), # The time drawn by the treatment as they stop (see post_art_dwell())
            hiv_death = ss.Transition(off_art, ss.dead, rate=self.death_rate),
            art_death = ss.Transition(on_art, ss.dead, rate=self.art_death_rate),
        )

        # Stopping treatment early, at the treatment's rate (sti.ART's rate_art_dropout); sti.ART.stop_art() takes people off it otherwise
        self.define_transitions(
            art_dropout = ss.Transition(on_art, self.art_discontinued, rate=self.dropout_rate),
        )

        # Results of the flows
        self.define_results(
            new_deaths = self.aids_death + self.post_art_death + self.hiv_death + self.art_death,
            new_infections_sex = 1*self.infection,
            new_infections_prenatal = 1*self.mtct_prenatal,
            new_infections_postnatal = 1*self.mtct_postnatal,
            new_infections_mtct = self.mtct_prenatal + self.mtct_postnatal,
        )

        # As code: the CD4 count, diagnosis, PrEP, and circumcision
        self.define_states(
            ss.FloatArr('ti_infected'), # Time of infection, as it was
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
        self.treatment_name = None # The module that treats HIV (e.g. sti.ART), which registers itself (see add_treatment())
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

            # CD4 count
            cd4_start = ss.normal(mean=800, std=50),
            cd4_latent = ss.normal(mean=500, std=50),

            # Infectiousness and susceptibility
            eff_condom = 0.9, # Read by the sexual networks (see sti.BaseNetwork.net_beta())
            rel_trans_acute = ss.normal(mean=5.3, std=0.5), # Bellan 2015: RH=5.3, EHM≈7.3
            rel_trans_falling = ss.normal(mean=8, std=0.5), # Increase in late infection (a CD4 count below 200)
            rel_sus_age = None, # Age/sex-stratified rel_sus multipliers: (age_lo, age_hi, sex, multiplier) tuples; sex is 'f', 'm', or None for both. E.g. [(15, 25, 'f', 1.7), (25, 50, 'f', 1.0), (15, 50, 'm', 1.0)]
            eff_circ = 0.6, # Program (medical) VMMC: reduction in HIV acquisition susceptibility
            eff_circ_traditional = 0.6, # Traditional (non-program) circumcision: reduction in HIV acquisition susceptibility

            # Initial conditions
            rel_init_prev = 1, # Scales init_prev_data
            init_diagnosed = ss.bernoulli(p=0),
            dist_ti_init_infected = ss.uniform(low=-10*12, high=-5), # When the infections at the start began, in steps

            # Treatment: sti.ART's parameters, which it takes from here if they are given to HIV (rate_art_dropout is the rate of the transition art_dropout, which is the treatment's)
            **{key: None for key in art_pars if key != 'rate_art_dropout'},
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

    # Treatment: the interface for the module that treats HIV (e.g. sti.ART)
    def add_treatment(self, treatment):
        """
        Register the module that treats HIV, which moves people into and out of the states of treatment
        (`on_effective_art`, `on_nonsuppressive_art`, `art_discontinued`). HIV calls it for: `update_cd4()`, the CD4
        counts of those it treats or has treated; `on_art_death_rate(uids)`, the annual rate of death of those on it;
        `dropout_rate(uids)`, the rate at which they stop it early; `post_art_dwell(uids)`, as they stop, the time until
        their CD4 count falls to zero; `trans_multiplier(uids)`, the multiplier of the infectiousness of those on it;
        and `refer(uids, ti)`, to schedule a start after a diagnosis
        """
        if self.treatment_name not in [None, treatment.name]:
            errormsg = f'HIV is treated by "{self.treatment_name}" already, so cannot be treated by "{treatment.name}" too: a sim has one sti.ART'
            raise ValueError(errormsg)
        self.treatment_name = treatment.name
        return

    @property
    def treatment(self):
        """ The module that treats HIV (e.g. sti.ART), if any """
        return self.sim.interventions[self.treatment_name] if self.treatment_name is not None else None

    @property
    def art(self):
        """ The module that holds the states of ART (sti.ART), if any """
        return self.treatment

    def refer(self, uids, ti=None):
        """ Refer people diagnosed to treatment, to start at `ti` (by default now); nothing if there is no treatment """
        if self.treatment is not None:
            self.treatment.refer(ss.uids(uids), self.ti if ti is None else ti)
        return

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

    def given_art_pars(self):
        """ The parameters of sti.ART given to HIV, which ART takes (see sti.ART.init_pre()) """
        placeholders = {key: self.pars[key] for key in art_pars if key != 'rate_art_dropout' and self.pars[key] is not None}
        return sc.mergedicts(placeholders, self.art_given)

    def death_rate(self, sim, uids):
        """ The rate of the declared transition hiv_death """
        return ss.peryear(self.off_art_death_rate(uids))

    def art_death_rate(self, sim, uids):
        """ The rate of the declared transition art_death: the treatment's (see sti.ART.on_art_death_rate()) """
        tx = self.treatment
        return ss.peryear(tx.on_art_death_rate(uids) if tx is not None else np.zeros(len(uids)))

    def dropout_rate(self, sim, uids):
        """ The rate of the declared transition art_dropout: the treatment's (see sti.ART.dropout_rate()) """
        tx = self.treatment
        return tx.dropout_rate(uids) if tx is not None else ss.peryear(np.zeros(len(uids)))

    def post_art_dwell(self, sim, uids):
        """ As people stop treatment (by sti.ART.stop_art() or art_dropout): the time until their CD4 count falls to zero, when they die of AIDS (post_art_death), from the treatment """
        return self.treatment.post_art_dwell(uids)

    def init_post(self):
        """ CD4 counts; and the infections at the start, which began before it """
        self.init_cd4()
        super().init_post()
        seeds = self.acute.uids
        self.backdate_seeds(seeds)
        self.set_infection(seeds)
        diagnosed = self.pars.init_diagnosed.filter(seeds)
        self.diagnosed[diagnosed] = True
        self.ti_diagnosed[diagnosed] = 0
        self.refer(diagnosed, ti=0) # Scheduled to start ART at once, so an ART intervention picks them up
        unused = list(self.given_art_pars()) # Taken by sti.ART when there is one (see its init_pre())
        if unused:
            errormsg = f'ART parameters were given to HIV ({sc.strjoin(unused)}), but there is no ART in the sim to use them: add sti.ART()'
            raise ValueError(errormsg)
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
        steps = lambda tr: eng.whole_steps(uids, eng.sample_dwell(tr, uids))
        to_latent = began + steps(trs.latency)
        to_falling = to_latent + steps(trs.progression)
        to_zero = to_falling + steps(trs.aids_death)
        self.ti_infected[uids] = began
        for state, start, end in [('acute', began, to_latent), ('latent', to_latent, to_falling), ('falling', to_falling, to_zero)]:
            here = (start <= now) & ((end > now) | (state == 'falling')) # Those whose course ended before the start are due to die in its first step
            if here.any():
                self.move(uids[here], self[state])
                self[f'ti_{state}'].raw[uids[here]] = start[here] # [PROPOSAL ST-E] ss.move(..., since=): when they entered it, and when they leave it (by now at the latest)
                self[f'ti_exit_{state}'].raw[uids[here]] = np.maximum(end[here], now)
        return

    # CD4 count
    def init_cd4(self):
        """ Set the CD4 counts of new agents """
        uids = ss.uids(np.isnan(self.cd4_start.raw).nonzero()[-1])
        self.cd4_start[uids] = self.pars.cd4_start.rvs(uids)
        self.cd4_nadir[uids] = sc.dcp(self.cd4_start[uids])
        return

    def progress(self, state, uids):
        """ The share of their time in a declared state that agents have spent, from 0 to 1 """
        return self.engine.clock_values(self.graph.states[state].progress, uids)

    def update_cd4(self):
        """ CD4 counts: falling from the start of infection to the latent level in acute infection, and from it to 1 in late infection; and those of the people in the states of treatment, from the treatment """
        acute = self.acute.uids
        cd4_start, cd4_latent = self.cd4_start[acute], self.cd4_latent[acute]
        self.cd4[acute] = cd4_start - (cd4_start - cd4_latent)*self.progress('acute', acute)
        latent = self.latent.uids
        self.cd4[latent] = self.cd4_latent[latent]
        falling = self.falling.uids
        self.cd4[falling] = np.maximum(1, self.cd4_latent[falling] - (self.cd4_latent[falling] - 1)*self.progress('falling', falling))
        if self.treatment is not None:
            self.treatment.update_cd4()
        untreated = (self.infected & ~self.on_art).uids
        self.cd4_nadir[untreated] = np.minimum(self.cd4_nadir[untreated], self.cd4[untreated])
        return

    def update_transmission(self):
        """
        Each agent's susceptibility (age, PrEP, circumcision), reset at the start of each step; connectors change it after this, and
        infectiousness, reset with it, is set by stage once the declared transitions have been taken (see `stage_infectiousness()`)
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

        return

    def stage_infectiousness(self):
        """ Infectiousness in acute and late infection, by the stage people are in after the step's declared transitions (a multiplier, as the connectors' are) """
        p = self.pars
        self.rel_trans[self.acute] *= p.rel_trans_acute.rvs(self.acute.uids)
        self.rel_trans[self.aids] *= p.rel_trans_falling.rvs(self.aids.uids)
        return

    def treat_infectiousness(self):
        """ The infectiousness of those on treatment, times the treatment's multiplier (see sti.ART.trans_multiplier()); applied just before transmission, after the interventions that start treatment in the step """
        uids = self.on_art.uids
        if len(uids) and self.treatment is not None:
            self.rel_trans[uids] *= self.treatment.trans_multiplier(uids)
        return

    def infectiousness(self, uids=None):
        """ How infectious each agent is, all told, relative to untreated infection (their rel_trans, which includes treatment's effect) """
        uids = self.sim.people.auids if uids is None else ss.uids(uids)
        return self.rel_trans[uids]

    def step_state(self):
        """ At the start of the step: CD4 counts, and each agent's susceptibility (the declared transitions are taken in step()) """
        self.init_cd4()
        self.update_cd4()
        self.update_transmission()
        return

    def step(self):
        """ The declared transitions; infectiousness by stage, and with ART's rise to full efficacy; transmission; then the CD4 counts of those infected (a child's from its mother) """
        eng = self.engine
        eng.step_state() # [PROPOSAL ST-F] The transitions not by contact: done here, so that infectiousness can follow them before transmission
        died = self.aids_death.uids | self.post_art_death.uids | self.hiv_death.uids | self.art_death.uids
        self.ti_dead[died] = self.ti
        self.stage_infectiousness()
        self.treat_infectiousness()
        eng.step() # Transmission
        eng.step_funcs()
        if self.callbacks:
            self.run_callbacks()
        children = self.mtct_prenatal.uids | self.mtct_postnatal.uids
        mothers, _ = self.engine.sources(children)
        known = ~np.isnan(mothers)
        self.cd4_start[children[known]] = self.cd4_start[ss.uids(mothers[known])]
        self.cd4_nadir[children] = self.cd4_start[children]
        infections = [self[name].uids for name,tr in self.graph.transitions.items() if 'acute' in [state.name for state in tr.dest]] # However they were infected (e.g. by a transition a subclass adds)
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

    def start_prep(self, uids, eff, dur, source_id, adh=1.0):
        """
        Start a PrEP course for uids not already on PrEP

        Args:
            uids: agents to start
            eff: base efficacy (0-1) of this course when fully adherent
            dur: course duration (an ss.Duration/Time value) before renewal is needed
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
            ss.Result('infected', dtype=int, kind='stock', label='Number infected'),
            ss.Result('on_art', dtype=int, kind='stock', label='On ART'),
            ss.Result('diagnosed', kind='stock', dtype=int, label='Diagnosed', auto_plot=False),
            ss.Result('cum_deaths', dtype=int, label='Cumulative deaths', auto_plot=False),
            ss.Result('new_diagnoses', dtype=int, label='Diagnoses'),
            ss.Result('cum_diagnoses', dtype=int, label='Cumulative diagnoses', auto_plot=False),
            ss.Result('new_agents_on_art', dtype=int, label='New treated', auto_plot=False),
            ss.Result('prevalence_15_49', kind='ratio', dtype=float, label='Prevalence 15-49', scale=False),
            ss.Result('prevalence_sw', kind='ratio', dtype=float, label='FSW prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_sw', dtype=int, label='New infections - FSW', auto_plot=False),
            ss.Result('new_infections_not_sw', dtype=int, label='New infections - Other F', auto_plot=False),
            ss.Result('prevalence_client', kind='ratio', dtype=float, label='Client prevalence', scale=False, auto_plot=False),
            ss.Result('new_infections_client', dtype=int, label='New infections - Clients', auto_plot=False),
            ss.Result('new_infections_not_client', dtype=int, label='New infections - Other M', auto_plot=False),
            ss.Result('p_on_art', kind='ratio', dtype=float, label='Proportion on ART', scale=False),
        ]
        if self.include_mtct:
            results += [
                ss.Result('n_on_art_pregnant', dtype=int, auto_plot=False),
                ss.Result('p_diagnosed_pregnant', kind='ratio', dtype=float, label='Proportion of HIV+ pregnant women diagnosed', scale=False, auto_plot=False),
            ]
        if 'structuredsexual' in self.sim.networks.keys():
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    results += [
                        ss.Result('prevalence_risk_group_' + str(risk_group) + '_' + sex, kind='ratio', scale=False, auto_plot=False),
                        ss.Result('new_infections_risk_group_' + str(risk_group) + '_' + sex, dtype=int, auto_plot=False),
                    ]
        self.define_results(*[res for res in results if res.name not in self.results and res.name not in self.graph.exprs and res.name not in self.graph.outcomes])
        return

    def record_stocks(self):
        super().record_stocks()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        t = self.tally()
        self.record_sti_stocks(t)
        infected, on_art, diagnosed = t.vals(self.infected), t.vals(self.on_art), t.vals(self.diagnosed)
        res['prevalence'][ti] = t.div(t.count(infected), t.count(t.vals(ppl.alive))) # For the whole population
        res['on_art'][ti] = t.count(on_art)
        res['diagnosed'][ti] = t.count(diagnosed)
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
            risk_group = np.asarray(nw.risk_group.values)
            for rg in range(nw.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    in_group = (risk_group == rg) & t.vals(ppl[sex])
                    res[f'prevalence_risk_group_{rg}_{sex}'][ti] = t.div(t.count(infected & in_group), t.count(in_group))
        return

    def record_flows(self):
        super().record_flows()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        t = self.tally()
        self.record_sti_flows(t)
        new_inf, on_art = t.vals(self.new_infected()), t.vals(self.on_art)
        res['new_diagnoses'][ti] = t.count(t.vals(self.ti_diagnosed == ti))
        if self.treatment is not None:
            res['new_agents_on_art'][ti] = t.count(t.vals(self.treatment.ti_art == ti) & on_art) # Gated on on_art, since ti_art may also hold a scheduled future start

        # By sex work and risk group
        if 'structuredsexual' in self.sim.networks.keys():
            nw = self.sim.networks.structuredsexual
            for group, key in [(t.vals(nw.fsw), 'sw'), (t.vals(nw.client), 'client')]:
                if t.count(group):
                    res[f'new_infections_{key}'][ti] = t.count(new_inf & group)
                    res[f'new_infections_not_{key}'][ti] = t.count(new_inf & ~group)
            risk_group = np.asarray(nw.risk_group.values)
            for rg in range(nw.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    in_group = (risk_group == rg) & t.vals(ppl[sex])
                    res[f'new_infections_risk_group_{rg}_{sex}'][ti] = t.count(new_inf & in_group)
        return

    def finalize_results(self):
        super().finalize_results()
        self.finalize_sti_results()
        res = self.results
        res['cum_deaths'][:] = np.cumsum(res['new_deaths'][:])
        res['cum_diagnoses'][:] = np.cumsum(res['new_diagnoses'][:])
        return



def _art_attr(name):
    def get(self):
        art = self.art
        if art is None:
            errormsg = f'hiv.{name} is a state or method of treatment, held by sti.ART, and there is no ART in this sim'
            raise AttributeError(errormsg)
        return getattr(art, name)
    get.__doc__ = f'The ART module\'s `{name}` (see sti.ART)'
    return property(get)

for _name in art_attrs:
    setattr(HIV, _name, _art_attr(_name))
