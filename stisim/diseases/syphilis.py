"""
Define syphilis disease module
"""

import numpy as np
import sciris as sc
from sciris import randround as rr # Since used frequently
import starsim as ss
import stisim as sti
from stisim.diseases.sti import STIResults, default_age_bins, default_sex_keys

__all__ = ['Syphilis', 'SyphilisPlaceholder']


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

    @property
    def symptomatic(self):
        """ Alias for the active state, so connectors share an interface with the full Syphilis module """
        return self.active

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


# Probabilities of the birth outcomes of a congenital infection, by the mother's stage when it happened:
# miscarriage, neonatal death, stillbirth, congenital syphilis, and a live birth without complications (which may be
# preterm or of low birth weight). Sources: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5973824/, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819963/
default_birth_outcomes = sc.objdict(
    mat_active = [0.00, 0.25, 0.25, 0.45, 0.05], # Mothers with primary or secondary infection
    early      = [0.00, 0.05, 0.10, 0.40, 0.45], # Mothers with early latent infection
    late       = [0.00, 0.00, 0.10, 0.10, 0.80], # Mothers with late latent infection
)
fetal_outcomes = ['miscarriage', 'nnd', 'stillborn', 'congenital', 'normal'] # In the order of the birth outcomes


class Syphilis(ss.Infection, STIResults):
    """
    Syphilis: exposed, primary, secondary, early and late latent, and tertiary stages, and congenital syphilis

    People infected are exposed (incubating) for `dur_exposed`, then have a primary chancre for `dur_primary`, then
    secondary syphilis for `dur_secondary`, then latent syphilis: early for `dur_early`, then late. A share of the
    latent (`p_reactivate`) relapse to secondary `time_to_reactivate` after the latent stage began; of the rest, a
    share (`p_tertiary`) progress to tertiary syphilis `time_to_tertiary` after it began, and of those a share
    (`p_death`) die `time_to_death` later. Infectiousness is by stage, and in the latent stages decays with a
    half-life of `rel_trans_latent_half_life` from the start of the latent stage.

    A fetus is infected in pregnancy at a rate `beta_m2c` (per month) while its mother is infectious, and its outcome
    (miscarriage, stillbirth, neonatal death, congenital syphilis, or none) depends on the mother's stage when it was
    infected (`birth_outcomes`). Fetal losses happen at the delivery (in the step before it, while the fetus is unborn:
    `ss.Pregnancy` steps before the diseases), and neonatal deaths and congenital syphilis at the birth.

    Serology (treponemal and non-treponemal test status) and whether symptoms are visible are each person's own,
    updated as the stages change.

    States: `susceptible`, `exposed`, `primary`, `secondary`, `early`, `late`, `tertiary`, the fetal states
    (`fetus_miscarriage`, `fetus_nnd`, `fetus_stillborn`, `fetus_congenital`, `fetus_normal`: infected in utero,
    by the outcome to come), and `congenital` (born with congenital syphilis); `infected` is any adult stage.
    `latent` (`early | late`), `infectious`, `symptomatic`, and the others below are conditions over them.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'syph'``)
        birth_outcomes (dict): the probabilities of the five birth outcomes for a mother with active (`mat_active`), `early` latent, and `late` latent syphilis, each a list or an `ss.choice(a=5, p=...)`
        init_latent_prev (float): the share of people who start in early latent syphilis
        age_range (list): age range of the main (adult) results (default [15, 50])
        age_bins (list): bin edges for age-stratified results; ``None`` skips them
        sex_keys (dict): sex strata for results; ``None`` skips sex-stratified results
        include_care (bool): whether to define the results of care
        kwargs (dict): parameter values (e.g. ``beta_m2f=0.2``), `init_prev`, or `label`
    """
    stages = ['exposed', 'primary', 'secondary', 'early', 'late', 'tertiary']
    fetal_states = [f'fetus_{o}' for o in fetal_outcomes]

    def __init__(self, pars=None, name='syph', birth_outcomes=None, init_latent_prev=None, age_range=None,
                 age_bins=default_age_bins, sex_keys=default_sex_keys, include_care=True, label=None, **kwargs):
        super().__init__(name=name, label=label)
        self.define_pars(
            dt = 'month',
            init_prev = 0.0,

            # Transmission: a probability per act on the sexual networks, by stage, and a rate per month from mother to fetus
            beta = 0.0, # On each sexual network, from the ones below (see validate_beta())
            beta_m2f = 0.1, # From a man to a woman
            rel_beta_f2m = 0.5, # From a woman to a man, relative to from a man to a woman
            beta_m2m = 0.0, # Between men
            beta_m2c = 0.0, # From mother to fetus in pregnancy, per month
            eff_condom = 0.0, # Read by the sexual networks (see sti.BaseNetwork.net_beta())
            rel_trans_primary = 1.0, # Infectiousness in the primary stage
            rel_trans_secondary = 1.0, # In the secondary stage, and at the start of the latent stage
            rel_trans_tertiary = 0.0, # In the tertiary stage
            rel_trans_latent_half_life = ss.years(1),

            # Natural history
            dur_exposed = ss.normal(ss.days(50), ss.days(10)), # Incubation period: https://pubmed.ncbi.nlm.nih.gov/9101629/
            dur_primary = ss.normal(ss.weeks(6), ss.weeks(1)), # https://pubmed.ncbi.nlm.nih.gov/9101629/
            dur_secondary = ss.lognorm_ex(ss.months(3.6), ss.months(1.5)), # https://pubmed.ncbi.nlm.nih.gov/9101629/
            dur_early = ss.uniform(ss.months(22), ss.months(24)), # Matches WHO Europe's definition (up to 2 years after infection)
            p_reactivate = 0.35, # Probability of relapsing from latent to secondary syphilis
            time_to_reactivate = ss.lognorm_ex(ss.years(1), ss.years(1)), # From the start of the latent stage
            p_tertiary = 0.35, # Probability of tertiary syphilis, for those who do not relapse: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917057/
            time_to_tertiary = ss.normal(ss.years(20), ss.years(2)), # From the start of the latent stage
            p_death = 0.05, # Probability of dying of tertiary syphilis
            time_to_death = ss.lognorm_ex(ss.years(5), ss.years(5)), # From tertiary syphilis

            # Visible symptoms, [female, male]: primary chancres are internal in women (less visible), external in men; ~60% go unnoticed
            p_symp_primary = [0.3, 0.5],
            p_symp_secondary = [0.85, 0.85], # A disseminated rash, visible in both sexes

            # Serology
            time_to_undetectable = ss.lognorm_ex(ss.years(15), ss.years(5)), # From late latent until the non-treponemal titre is below detection
            dur_to_trep = ss.normal(ss.weeks(3), ss.weeks(1)), # Window period: from infection until treponemal antibodies are detected
            p_trep_persists = ss.bernoulli(p=0.80), # The other 20% sero-revert eventually
            time_to_trep_revert = ss.normal(ss.years(10), ss.years(5)),
            anc_detection = 0.8,
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            *[ss.BoolState(name) for name in self.stages + self.fetal_states + ['congenital']],
            ss.BoolState('ever_exposed'), # Ever infected
            ss.IntArr('n_infections', default=0), # Infections in a lifetime
            ss.BoolState('trep'), # Treponemal test positive: 80% stay so for life, 20% sero-revert eventually
            ss.BoolState('nontrep'), # Non-treponemal (RPR/VDRL) test positive
            ss.BoolState('chancre_visible'), # Whether the primary chancre is at a visible site
            ss.BoolState('rash_visible'), # Whether the secondary rash is noticed
            ss.TimeArr('ti_primary'), # When each stage starts, drawn when the one before it does (NaN for never)
            ss.TimeArr('ti_secondary'), # From primary, or a relapse from latent
            ss.TimeArr('ti_latent'), # When the latent stage began, for the decay of infectiousness
            ss.TimeArr('ti_late'),
            ss.TimeArr('ti_tertiary'),
            ss.TimeArr('ti_dead'),
            ss.TimeArr('ti_outcome'), # When the outcome of a fetus infected in utero happens
            ss.TimeArr('ti_nontrep_end'), # When the non-treponemal titre falls below detection
            ss.TimeArr('ti_trep_start'), # When treponemal antibodies appear; NaN if treated in the window period
            ss.TimeArr('ti_trep_end'), # When the treponemal test reverts (NaN for those for whom it persists)
        )
        outcomes = sc.mergedicts(default_birth_outcomes, birth_outcomes)
        self.birth_outcomes = {stage:np.cumsum(val.pars.p if isinstance(val, ss.choice) else val) for stage,val in outcomes.items()} # Cumulative, also from Starsim 3's ss.choice
        self.init_latent_prev = init_latent_prev
        self.u = sc.objdict({key:ss.random(name=f'u_{key}') for key in ['chancre', 'rash', 'reactivate', 'tertiary', 'death', 'mtct', 'outcome', 'init_latent']}) # A uniform for each chance, so each person's is their own
        self.set_strata(age_range if age_range is not None else [15, 50], age_bins, sex_keys)
        self.include_care = include_care
        self.store_sw = False # Whether to store results for sex workers and their clients
        self.store_risk_groups = False # Whether to store results by risk group
        return

    @property
    def latent(self):
        """ Latent: early or late """
        return self.early | self.late

    @property
    def infectious(self):
        """ Infectious: primary, secondary, or latent (and tertiary, whose infectiousness is `rel_trans_tertiary`, 0 by default) """
        return self.symptomatic | self.latent | self.tertiary

    @property
    def treatable(self):
        """ Who a treatment can cure: everyone infected """
        return self.infected

    @property
    def naive(self):
        """ Never exposed """
        return ~self.ever_exposed

    @property
    def sus_not_naive(self):
        """ Susceptible but with syphilis antibodies, which persist after treatment """
        return self.susceptible & self.ever_exposed

    @property
    def mat_active(self):
        """ Definition of infection used for maternal outcomes, includes exposed + primary + secondary stages """
        return self.exposed | self.primary | self.secondary

    @property
    def sexually_transmissible(self):
        """ Stages during which syphilis is sexually transmissible: primary,
        secondary, and early latent. Late latent and tertiary are not
        normally sexually transmissible (late latent can still transmit
        congenitally at low levels). Matches WHO 'early infectious syphilis'. """
        return self.primary | self.secondary | self.early

    @property
    def symptomatic(self):
        """ Symptomatic stages of syphilis: primary or secondary. Distinct from
        ``visibly_symptomatic``, which requires that the symptom actually be
        observable (chancre at a visible site, or rash). """
        return self.primary | self.secondary

    @property
    def ulcerative(self):
        """ Has a visible genital ulcer (chancre at a visible site during primary stage) """
        return self.chancre_visible & self.primary

    @property
    def visibly_symptomatic(self):
        """ Has clinically detectable symptoms: visible ulcer or visible rash.
        Used by care-seeking pathways. """
        return self.ulcerative | (self.rash_visible & self.secondary)

    @property
    def fetal(self):
        """ Infected in utero, whatever the outcome to come """
        return self.fetus_miscarriage | self.fetus_nnd | self.fetus_stillborn | self.fetus_congenital | self.fetus_normal

    def validate_beta(self):
        """ The betas of the sexual networks, from the parameters; the maternal network does not transmit, since a fetus is infected at a rate of its own (see step()) """
        p = self.pars
        self.pars.beta = sti.sexual_beta(p.beta_m2f, p.beta_m2f*p.rel_beta_f2m, p.beta_m2m)
        return super().validate_beta()

    def by_sex(self, pair, uids):
        """ The value of a [female, male] pair for each of these people """
        return np.where(self.sim.people.female.raw[uids], pair[0], pair[1])

    def chance(self, key, p, uids):
        """ Whether each of some agents has an event of probability p (a number, or one for each) """
        return self.u[key].rvs(uids) < p

    def new_infected(self):
        """ Who was infected (not in utero) in this step, including those seeded at the start, as a boolean array """
        return self.ti_infected == self.ti

    def init_post(self):
        """ The cases at the start (see `ss.Infection.init_post()`), and those who start in early latent syphilis """
        super().init_post()
        if self.init_latent_prev:
            uids = self.susceptible.uids
            uids = uids[self.chance('init_latent', self.init_latent_prev, uids)]
            self.set_prognoses(uids)
            self.exposed[uids] = False
            self.nontrep[uids] = True # As after the secondary stage
            self.set_latent(uids)
        return

    def set_prognoses(self, uids, sources=None):
        """ Infect some adults: their incubation and primary stage, serology to come, and whether the chancre will be visible """
        super().set_prognoses(uids, sources)
        ti = self.ti
        p = self.pars
        self.infected[uids] = True
        self.exposed[uids] = True
        self.ti_infected[uids] = ti if self.sim.initialized else 0 # The infections seeded at the start are of its first point (see SEIS.set_prognoses())
        self.ever_exposed[uids] = True
        self.n_infections[uids] += 1
        start = self.ti + 1 if self.sim.initialized else self.ti # A course runs from when the infection happens: at the end of the step, in transmission (as a declared disease dates it), or as ss.Infection.init_post() dates those seeded at the start
        self.ti_primary[uids] = start + p.dur_exposed.rvs(uids)
        self.ti_secondary[uids] = self.ti_primary[uids] + p.dur_primary.rvs(uids)

        # Treponemal seroconversion after the window period; 80% stay positive for life, the others sero-revert later
        self.ti_trep_start[uids] = start + p.dur_to_trep.rvs(uids)
        reverters = uids[~p.p_trep_persists.rvs(uids)]
        self.ti_trep_end[reverters] = self.ti_trep_start[reverters] + p.time_to_trep_revert.rvs(reverters)
        self.chancre_visible[uids] = self.chance('chancre', self.by_sex(p.p_symp_primary, uids), uids)
        return

    def set_secondary(self, uids):
        """ Entering the secondary stage (from primary, or a relapse): its end, the non-treponemal titre rises, and the rash may be noticed """
        p = self.pars
        for state in ['primary', 'early', 'late']:
            self[state][uids] = False
        self.secondary[uids] = True
        self.ti_latent[uids] = self.ti + p.dur_secondary.rvs(uids)
        self.nontrep[uids] = True
        self.rash_visible[uids] = self.chance('rash', self.by_sex(p.p_symp_secondary, uids), uids)
        return

    def set_latent(self, uids):
        """ Entering the latent stage: the rash clears and infectiousness begins to decay; the time to late latent, and a relapse, or tertiary syphilis and death, to come """
        ti = self.ti
        p = self.pars
        self.secondary[uids] = False
        self.early[uids] = True
        self.rash_visible[uids] = False
        self.ti_latent[uids] = ti
        self.ti_late[uids] = ti + p.dur_early.rvs(uids)
        relapse = self.chance('reactivate', p.p_reactivate, uids)
        self.ti_secondary[uids] = np.nan
        self.ti_secondary[uids[relapse]] = ti + p.time_to_reactivate.rvs(uids[relapse])
        rest = uids[~relapse]
        tert = rest[self.chance('tertiary', p.p_tertiary, rest)]
        self.ti_tertiary[tert] = ti + p.time_to_tertiary.rvs(tert)
        die = tert[self.chance('death', p.p_death, tert)]
        self.ti_dead[die] = self.ti_tertiary[die] + p.time_to_death.rvs(die)
        return

    def step_state(self):
        """ Take people through the stages whose time has come; then update serology, fetal outcomes, and the infectiousness of each stage """
        ti = self.ti
        p = self.pars
        due = lambda state, when: (state & when.due()).uids

        self.infect_fetuses()
        uids = due(self.exposed, self.ti_primary)
        self.exposed[uids] = False
        self.primary[uids] = True
        self.set_secondary(due(self.primary, self.ti_secondary) | due(self.latent, self.ti_secondary)) # From primary, or a relapse
        self.set_latent(due(self.secondary, self.ti_latent))
        uids = due(self.early, self.ti_late)
        self.early[uids] = False
        self.late[uids] = True
        self.ti_nontrep_end[uids] = ti + p.time_to_undetectable.rvs(uids)
        uids = due(self.latent, self.ti_tertiary)
        self.early[uids] = self.late[uids] = False
        self.tertiary[uids] = True
        dead = due(self.tertiary, self.ti_dead)
        self.sim.people.request_death(dead)
        self.results.new_deaths[ti] = len(dead)

        # Fetal outcomes: losses (in the step before the delivery), and neonatal deaths, congenital syphilis, and none at the birth
        res = self.results
        for outcome, key in zip(fetal_outcomes, ['new_miscarriages', 'new_nnds', 'new_stillborns', 'new_congenital', None]):
            uids = due(self[f'fetus_{outcome}'], self.ti_outcome)
            self[f'fetus_{outcome}'][uids] = False
            self.ti_outcome[uids] = np.nan
            if outcome == 'congenital':
                self.congenital[uids] = True
            elif outcome == 'normal':
                self.susceptible[uids] = True
            else:
                self.sim.people.request_death(uids)
            if key:
                res[key][ti] = len(uids)
        res.new_congenital_deaths[ti] = res.new_miscarriages[ti] + res.new_nnds[ti] + res.new_stillborns[ti]

        # Serology over time
        self.nontrep[due(self.nontrep, self.ti_nontrep_end)] = False
        self.trep[due(~self.trep, self.ti_trep_start)] = True
        self.trep[due(self.trep, self.ti_trep_end)] = False

        # Susceptibility and infectiousness of each person (which connectors change after this), with the decay in the latent stage
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1
        self.rel_trans[self.primary] = p.rel_trans_primary
        self.rel_trans[self.secondary] = p.rel_trans_secondary
        self.rel_trans[self.tertiary] = p.rel_trans_tertiary
        latent = self.latent.uids
        self.rel_trans[latent] = p.rel_trans_secondary*0.5**((ti - self.ti_latent[latent])/p.rel_trans_latent_half_life.to_dt())
        return

    def step(self):
        """ Sexual transmission (see `ss.Infection.step()`) """
        new, sources, networks = super().step()
        self.sexual_sources = sources # All of its transmission over networks is sexual
        self.results.new_infections_sex[self.ti] = len(new)
        return new, sources, networks

    def infect_fetuses(self):
        """ At the start of the step (as a rate, with the mothers' stages and the pregnancies then): fetuses of infectious mothers are infected at a rate, with an outcome by the mother's stage """
        res = self.results
        ti = self.ti
        res.new_infections_prenatal[ti] = res.new_infections_mtct[ti] = 0
        nets = [net for net in self.sim.networks.values() if isinstance(net, ss.PrenatalNet)]
        if not nets or not self.pars.beta_m2c or 'pregnancy' not in self.sim.demographics:
            return
        ppl = self.sim.people
        mothers = np.concatenate([np.asarray(net.edges.p1) for net in nets])
        fetuses = np.concatenate([np.asarray(net.edges.p2) for net in nets])
        at_risk = (ppl.age.raw[fetuses] < 0) & self.susceptible.raw[fetuses] & (self.symptomatic | self.latent).raw[mothers]
        mothers, fetuses = mothers[at_risk], ss.uids(fetuses[at_risk])
        hit = self.chance('mtct', -np.expm1(-self.pars.beta_m2c*self.t.dt_year*12), fetuses)
        mothers, fetuses = mothers[hit], fetuses[hit]
        if len(fetuses):
            u = self.u.outcome.rvs(fetuses)
            outcome = np.zeros(len(fetuses), dtype=int)
            for stage, cum in self.birth_outcomes.items():
                inside = (self.mat_active if stage == 'mat_active' else self[stage]).raw[mothers]
                outcome[inside] = np.minimum(np.searchsorted(cum, u[inside], side='right'), len(cum) - 1)
            delivery = np.ceil(self.sim.demographics.pregnancy.ti_exit_P.raw[mothers]) # The step of the delivery
            loss = np.isin(outcome, [0, 2])
            self.ti_outcome[fetuses] = np.where(loss, np.maximum(delivery - 1, ti), delivery) # Losses while the fetus is unborn, the rest once born
            for i, name in enumerate(self.fetal_states):
                self[name][fetuses[outcome == i]] = True
            self.susceptible[fetuses] = False
            res.new_infections_prenatal[ti] = res.new_infections_mtct[ti] = len(fetuses)
        return

    def clear_infection(self, uids):
        """ Cure people (e.g. treated), including congenital syphilis and infections in utero: they are susceptible again, and their serology stays as it is """
        for name in self.stages + self.fetal_states + ['infected', 'congenital']:
            self[name][uids] = False
        self.susceptible[uids] = True
        for name in ['primary', 'secondary', 'latent', 'late', 'tertiary', 'dead', 'outcome']:
            self[f'ti_{name}'][uids] = np.nan
        return

    def cure_fetuses(self, uids):
        """ Fetuses infected in utero (e.g. their mothers were treated) are born without complications """
        for name in self.fetal_states:
            self[name][uids] = name == 'fetus_normal'
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        results = [
            ss.Result('trep_prevalence', dtype=float, scale=False, label="Treponemal-test positive (ever-exposed proxy)"),
            ss.Result('pregnant_prevalence', dtype=float, scale=False, label="Pregnant prevalence", auto_plot=False),
            ss.Result('detected_pregnant_prevalence', dtype=float, scale=False, label="ANC prevalence", auto_plot=False),
            ss.Result('delivery_prevalence', dtype=float, scale=False, label="Delivery prevalence", auto_plot=False),
            ss.Result('sexually_transmissible_prevalence', dtype=float, scale=False, label="Sexually transmissible (primary + secondary + early latent)"),
            ss.Result('symptomatic_prevalence', dtype=float, scale=False, label="Symptomatic stage prevalence (primary + secondary)"),
            ss.Result('primary_prevalence', dtype=float, scale=False, label="Primary stage prevalence"),
            ss.Result('nontrep_prevalence', dtype=float, scale=False, label="Non-treponemal (RPR) positive"),
            ss.Result('n_active', dtype=int, label="Number in the primary or secondary stage", auto_plot=False),
            ss.Result('new_deaths', dtype=int, label="Deaths of tertiary syphilis", auto_plot=False),
            ss.Result('new_miscarriages', dtype=int, label="Congenital miscarriages", auto_plot=False),
            ss.Result('new_stillborns', dtype=int, label="Congenital stillbirths", auto_plot=False),
            ss.Result('new_nnds', dtype=int, label="Congenital neonatal deaths", auto_plot=False),
            ss.Result('new_congenital', dtype=int, label="Congenital cases", auto_plot=False),
            ss.Result('new_congenital_deaths', dtype=int, label="Congenital deaths", auto_plot=False),
            ss.Result('cum_congenital', dtype=int, label="Cumulative congenital cases", auto_plot=False),
            ss.Result('cum_congenital_deaths', dtype=int, label="Cumulative congenital deaths", auto_plot=False),
            ss.Result('new_reinfections', dtype=int, scale=True, label="Reinfections"),

            # Fetal testing and treatment, which might be assembled from several interventions
            ss.Result('new_fetus_treated_success', dtype=int, label="Fetal treatment success", auto_plot=False),
            ss.Result('new_fetus_treated_unnecessary', dtype=int, label="Fetal overtreatment", auto_plot=False),
            ss.Result('new_fetus_treated_failure', dtype=int, label="Fetal treatment failure", auto_plot=False),
            ss.Result('new_treated_unnecessary_pregnant', dtype=int, label="Overtreatment pregnant", auto_plot=False),
        ]

        # Most results are stored by age and sex
        for skk, skl, attr in self.sex_strata():
            if skk != '':
                results += [
                    ss.Result(f'trep_prevalence{skk}', scale=False, label=f"Treponemal-test positive{skl}", auto_plot=False),
                    ss.Result(f'sexually_transmissible_prevalence{skk}', scale=False, label=f"Sexually transmissible prevalence{skl}", auto_plot=False),
                    ss.Result(f'symptomatic_prevalence{skk}', scale=False, label=f"Symptomatic prevalence{skl}", auto_plot=False),
                    ss.Result(f'primary_prevalence{skk}', scale=False, label=f"Primary stage prevalence{skl}", auto_plot=False),
                    ss.Result(f'nontrep_prevalence{skk}', scale=False, label=f"Non-treponemal positive{skl}", auto_plot=False),
                ]
            for ab1, ab2 in self.age_strata():
                ask = f'{skk}_{ab1}_{ab2}'
                asl = f' ({skl}, {ab1}-{ab2})'
                results += [
                    ss.Result(f'sexually_transmissible_prevalence{ask}', scale=False, label=f"Sexually transmissible prevalence{asl}", auto_plot=False),
                    ss.Result(f'symptomatic_prevalence{ask}', scale=False, label=f"Symptomatic prevalence{asl}", auto_plot=False),
                    ss.Result(f'primary_prevalence{ask}', scale=False, label=f"Primary stage prevalence{asl}", auto_plot=False),
                    ss.Result(f'trep_prevalence{ask}', scale=False, label=f"Treponemal positive{asl}", auto_plot=False),
                    ss.Result(f'nontrep_prevalence{ask}', scale=False, label=f"Non-treponemal positive{asl}", auto_plot=False),
                ]

        # Add FSW and clients to results:
        if self.store_sw:
            results += [
                ss.Result('prevalence_sw', dtype=float, scale=False, label="Prevalence - FSW", auto_plot=False),
                ss.Result('new_infections_sw', dtype=int, label="Infections - FSW", auto_plot=False),
                ss.Result('new_infections_not_sw', dtype=int, label="Infections - other F", auto_plot=False),
                ss.Result('prevalence_client', dtype=float, scale=False, label="Prevalence - clients", auto_plot=False),
                ss.Result('new_infections_client', dtype=int, label="Infections - clients", auto_plot=False),
                ss.Result('new_infections_not_client', dtype=int, label="Infections - other M", auto_plot=False),
            ]

        # Add risk groups to results
        if self.store_risk_groups:
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
        new_inf = t.vals(self.new_infected())
        res.n_active[ti] = t.count(t.vals(self.symptomatic))

        # Prevalence among sexually active adults
        age = np.asarray(ppl.age.values)
        active = (age >= self.age_range[0]) & (age < self.age_range[1]) & t.vals(self.sim.networks.structuredsexual.active(ppl))
        stages = dict(prevalence=self.infected, trep_prevalence=self.trep, nontrep_prevalence=self.nontrep, sexually_transmissible_prevalence=self.sexually_transmissible,
                      symptomatic_prevalence=self.symptomatic, primary_prevalence=self.primary)
        stages = {key:t.vals(val) for key,val in stages.items()}
        for key, val in stages.items():
            res[key][ti] = t.div(t.count(val & active), t.count(active))

        # Pregnant women, if there is pregnancy
        if 'pregnancy' in self.sim.demographics.keys():
            infected, preg = stages['prevalence'], t.vals(ppl.pregnancy.P)
            res['pregnant_prevalence'][ti] = preg_prev = t.div(t.count(infected & preg), t.count(preg))
            res['detected_pregnant_prevalence'][ti] = preg_prev * self.pars.anc_detection
            delivered = t.vals(ppl.pregnancy.ti_delivered == ti)
            res['delivery_prevalence'][ti] = t.div(t.count(infected & delivered), t.count(delivered))
        res['new_reinfections'][ti] = t.count(new_inf & (np.asarray(self.n_infections.values) > 1))

        # Add FSW and clients to results:
        nw = self.sim.networks.structuredsexual
        if self.store_sw:
            for key, group in [('sw', t.vals(nw.fsw)), ('client', t.vals(nw.client))]:
                res[f'prevalence_{key}'][ti] = t.div(t.count(stages['prevalence'] & group), t.count(group))
                res[f'new_infections_{key}'][ti] = t.count(new_inf & group)
                res[f'new_infections_not_{key}'][ti] = t.count(new_inf & ~group)

        # Add risk groups
        if self.store_risk_groups:
            rg = np.asarray(nw.risk_group.values)
            for risk_group in range(nw.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    group = (rg == risk_group) & t.vals(ppl[sex])
                    if (new_inf & group).any():
                        res[f'prevalence_risk_group_{risk_group}_{sex}'][ti] = t.div(t.count(stages['prevalence'] & group), t.count(group))
                        res[f'new_infections_risk_group_{risk_group}_{sex}'][ti] = t.count(new_inf & group)

        # Results by age and sex
        for skk, skl, attr in self.sex_strata():
            sex = t.vals(ppl[attr])
            if skk != '':
                for key, val in stages.items():
                    res[f'{key}{skk}'][ti] = t.div(t.count(val & active & sex), t.count(active & sex)) # Age-binned prevalences (below): all people of that sex in each age bin
            if self.age_bins is not None:
                h = t.hist(sex)
                self.set_age_results(skk, **{key:t.hdiv(t.hist(stages[key] & sex), h) for key in ['sexually_transmissible_prevalence', 'symptomatic_prevalence', 'primary_prevalence', 'trep_prevalence', 'nontrep_prevalence']})
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        res['cum_congenital'][:] = np.cumsum(res['new_congenital'])
        res['cum_congenital_deaths'][:] = np.cumsum(res['new_congenital_deaths'])
        return
