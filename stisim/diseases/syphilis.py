"""
Define syphilis disease module
"""

import numpy as np
import sciris as sc
from sciris import randround as rr # Since used frequently
import starsim as ss
import stisim as sti
from stisim.diseases.sti import STIResults, default_age_bins, default_sex_keys
from stisim.utils import count, div, cond_prob

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


class Syphilis(ss.Disease, STIResults):
    """
    Syphilis: exposed, primary, secondary, early and late latent, and tertiary stages, and congenital syphilis, declared

    People infected are exposed (incubating) for `dur_exposed`, then have a primary chancre for `dur_primary`, then
    secondary syphilis for `dur_secondary`, then latent syphilis: early for `dur_early`, then late. A share of the
    latent (`p_reactivate`) relapse to secondary `time_to_reactivate` after the latent stage began; of the rest, a
    share (`p_tertiary`) progress to tertiary syphilis `time_to_tertiary` after it began, and of those a share
    (`p_death`) die `time_to_death` later. Infectiousness is by stage, and in the latent stages decays with a
    half-life of `rel_trans_latent_half_life` from the start of the latent stage (as agents).

    A fetus is infected in pregnancy at a rate `beta_m2c` (per month) while its mother is infectious, and its outcome
    at delivery (miscarriage, stillbirth, neonatal death, congenital syphilis, or none) depends on the mother's stage
    when it was infected (`birth_outcomes`). Fetal losses happen at the delivery (`ss.Pregnancy`'s `delivery`), and
    neonatal deaths and congenital syphilis at the birth.

    Serology (treponemal and non-treponemal test status) and whether symptoms are visible are not stages: they are
    arrays of each person's, updated at the events of the stages.

    States: `susceptible`, `exposed`, `primary`, `secondary`, `early`, `late`, `tertiary`, the fetal states
    (`fetus_miscarriage`, `fetus_nnd`, `fetus_stillborn`, `fetus_congenital`, `fetus_normal`: infected in utero,
    by the outcome to come), and `congenital` (born with congenital syphilis). `latent` (`early | late`),
    `infected`, `infectious`, `symptomatic`, and the others below are conditions over them.

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
    def __init__(self, pars=None, name='syph', birth_outcomes=None, init_latent_prev=None, age_range=None,
                 age_bins=default_age_bins, sex_keys=default_sex_keys, include_care=True, **kwargs):

        # Transmission: a probability per act on the sexual networks, by stage, and a rate per month from mother to fetus
        beta_m2f = ss.par(0.1, doc='Probability of transmission per act from a man to a woman')
        rel_beta_f2m = ss.par(0.5, doc='Transmission from a woman to a man, relative to from a man to a woman')
        beta_m2m = ss.par(0.0, doc='Probability of transmission per act between men')
        beta_m2c = ss.par(0.0, doc='Rate of transmission from mother to fetus in pregnancy, per month')
        rel_trans_primary = ss.par(1.0, doc='Infectiousness in the primary stage')
        rel_trans_secondary = ss.par(1.0, doc='Infectiousness in the secondary stage, and at the start of the latent stage')
        rel_trans_tertiary = ss.par(0.0, doc='Infectiousness in the tertiary stage')

        # Natural history
        dur_exposed = ss.par(ss.normal(ss.days(50), ss.days(10)), doc='Incubation period') # https://pubmed.ncbi.nlm.nih.gov/9101629/
        dur_primary = ss.par(ss.normal(ss.weeks(6), ss.weeks(1)), doc='Duration of the primary stage') # https://pubmed.ncbi.nlm.nih.gov/9101629/
        dur_secondary = ss.par(ss.lognorm_ex(ss.months(3.6), ss.months(1.5)), doc='Duration of the secondary stage') # https://pubmed.ncbi.nlm.nih.gov/9101629/
        dur_early = ss.par(ss.uniform(ss.months(22), ss.months(24)), doc='Duration of the early latent stage') # Matches WHO Europe's definition (up to 2 years after infection)
        p_reactivate = ss.par(0.35, doc='Probability of relapsing from latent to secondary syphilis')
        time_to_reactivate = ss.par(ss.lognorm_ex(ss.years(1), ss.years(1)), name='time_to_reactivate', doc='Time from the start of the latent stage to relapse')
        p_tertiary = ss.par(0.35, doc='Probability of progressing to tertiary syphilis, for those who do not relapse') # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917057/
        time_to_tertiary = ss.par(ss.normal(ss.years(20), ss.years(2)), name='time_to_tertiary', doc='Time from the start of the latent stage to tertiary syphilis')
        p_death = ss.par(0.05, doc='Probability of dying of tertiary syphilis')
        time_to_death = ss.par(ss.lognorm_ex(ss.years(5), ss.years(5)), doc='Time from tertiary syphilis to death')

        # States
        susceptible, exposed, primary, secondary, early, late, tertiary = ss.states('susceptible exposed primary secondary early late tertiary')
        fetal = ss.states('fetus_miscarriage fetus_nnd fetus_stillborn fetus_congenital fetus_normal') # Infected in utero, by the outcome to come
        fetus_miscarriage, fetus_nnd, fetus_stillborn, fetus_congenital, fetus_normal = fetal
        congenital = ss.State('congenital')
        latent = early + late

        # Transmission
        sexual = sti.sexual_beta(beta_m2f, beta_m2f*rel_beta_f2m, beta_m2m)
        infection = susceptible >> exposed @ sexual*(rel_trans_primary*primary + rel_trans_secondary*(secondary + latent) + rel_trans_tertiary*tertiary)

        # Natural history
        primary_onset = exposed >> primary @ dur_exposed
        secondary_onset = primary >> secondary @ dur_primary
        latency = secondary >> early @ dur_secondary
        to_late = early >> late @ dur_early
        reactivation = ss.Transition(latent, [secondary, ss.stay], p=[p_reactivate, None], on=latency, dur=time_to_reactivate) # Timed from the start of the latent stage, early or late
        tertiary_onset = ss.Transition(latent, [tertiary, ss.stay], p=[p_tertiary, None], on=latency, dur=time_to_tertiary)
        death = ss.Transition(tertiary, [ss.null, ss.stay], p=[p_death, None], dur=time_to_death)

        # Congenital syphilis: the outcome is decided by the mother's stage when the fetus is infected
        mat_active = ss.mother(exposed | primary | secondary) # As Starsim 3 classed the mothers
        mat_early = ss.mother(early)
        mat_late = ss.mother(late)
        outcomes = sc.mergedicts(default_birth_outcomes, birth_outcomes)
        probs = {stage:np.asarray(val.pars.p if isinstance(val, ss.choice) else val) for stage,val in outcomes.items()} # The probabilities, also from Starsim 3's ss.choice
        p_outcomes = [probs['mat_active'][i]*mat_active + probs['early'][i]*mat_early + probs['late'][i]*mat_late for i in range(4)] + [None] # A live birth without complications is the rest
        mtct = ss.Transition(susceptible, fetal, rate=ss.permonth(1)*beta_m2c*ss.mother(primary | secondary | latent), p=p_outcomes, acquisition=True)
        delivery = ss.state('pregnancy.delivery') # The events of ss.Pregnancy, which is made separately
        birth = ss.state('pregnancy.birth')
        miscarriage = ss.Transition(fetus_miscarriage, ss.null, on=delivery) # A fetal loss, at the delivery
        stillbirth = ss.Transition(fetus_stillborn, ss.null, on=delivery)
        nnd = ss.Transition(fetus_nnd, ss.null, on=birth) # A neonatal death, once born
        born_congenital = ss.Transition(fetus_congenital, congenital, on=birth)
        born_normal = ss.Transition(fetus_normal, susceptible, on=birth)

        # The parameters that the transitions do not read, which its own code does (below); and, as the graph does not yet register
        # them [PROPOSAL ST-D], the times after an event (on=), named above, whose values are set once it is built
        values = sc.mergedicts(dict(init_prev=0.0), pars, kwargs)
        code_pars = dict(
            eff_condom = 0.0, # Read by the sexual networks (see sti.BaseNetwork.net_beta())
            rel_trans_latent_half_life = ss.years(1),
            p_symp_primary = [0.3, 0.5], # [female, male]: primary chancres are internal in women (less visible), external in men; ~60% go unnoticed
            p_symp_primary_dist = ss.bernoulli(p=0), # Set for each person
            p_symp_secondary = [0.85, 0.85], # A disseminated rash, visible in both sexes
            p_symp_secondary_dist = ss.bernoulli(p=0),
            time_to_undetectable = ss.lognorm_ex(ss.years(15), ss.years(5)), # From late latent until the non-treponemal titre is below detection
            dur_to_trep = ss.normal(ss.weeks(3), ss.weeks(1)), # Window period: from infection until treponemal antibodies are detected
            p_trep_persists = ss.bernoulli(p=0.80), # The other 20% sero-revert eventually
            time_to_trep_revert = ss.normal(ss.years(10), ss.years(5)),
            anc_detection = 0.8,
        )
        delays = dict(time_to_reactivate=time_to_reactivate, time_to_tertiary=time_to_tertiary)
        given = {key:values.pop(key) for key in list(code_pars) + list(delays) if key in values}
        if init_latent_prev:
            values['init'] = sc.mergedicts(values.get('init'), dict(early=init_latent_prev))

        # [PROPOSAL ST-A] Declared here rather than in the class body, since on= a sibling transition does not yet work there
        super().__init__(
            infection=infection, primary_onset=primary_onset, secondary_onset=secondary_onset, latency=latency, to_late=to_late,
            reactivation=reactivation, tertiary_onset=tertiary_onset, death=death,
            mtct=mtct, miscarriage=miscarriage, stillbirth=stillbirth, nnd=nnd, born_congenital=born_congenital, born_normal=born_normal,
            n_active=primary + secondary, new_infections_sex=1*infection, new_infections_prenatal=1*mtct, new_infections_mtct=1*mtct, # Results of the states and flows
            new_nnds=1*nnd, new_stillborns=1*stillbirth, new_miscarriages=1*miscarriage, new_congenital=1*born_congenital,
            new_congenital_deaths=nnd + stillbirth + miscarriage, new_deaths=1*death,
            name=name, **sc.mergedicts(dict(dt='month'), values), # Monthly unless dt is given
        )
        self.define_pars(**code_pars, **delays)
        self.pars.update(given)

        # Not stages: serology and visible symptoms, each person's own, updated at the events of the stages (see step_state())
        self.define_states(
            ss.BoolState('ever_exposed'), # Ever infected
            ss.IntArr('n_infections', default=0), # Infections in a lifetime
            ss.BoolState('trep'), # Treponemal test positive: 80% stay so for life, 20% sero-revert eventually
            ss.BoolState('nontrep'), # Non-treponemal (RPR/VDRL) test positive
            ss.BoolState('chancre_visible'), # Whether the primary chancre is at a visible site
            ss.BoolState('rash_visible'), # Whether the secondary rash is noticed
            ss.FloatArr('ti_latent'), # When the latent stage began, for the decay of infectiousness
            ss.FloatArr('ti_nontrep_end'), # When the non-treponemal titre falls below detection
            ss.FloatArr('ti_trep_start'), # When treponemal antibodies appear; NaN if treated in the window period
            ss.FloatArr('ti_trep_end'), # When the treponemal test reverts (NaN for those for whom it persists)
        )

        # Results
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
    def infected(self):
        """ Infected (not in utero) """
        return self.exposed | self.primary | self.secondary | self.latent | self.tertiary

    @property
    def infectious(self):
        """ Infectious: primary, secondary, or latent """
        return self.symptomatic | self.latent

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

    def by_sex(self, pair, uids):
        """ The value of a [female, male] pair for each of these people """
        return np.where(self.sim.people.female[uids], pair[0], pair[1])

    def new_infected(self):
        """ Who was infected (not in utero) in this step: by transmission, or seeded at the start """
        new = self.infection.uids
        if self.ti == 0:
            new = new | self.infected.uids
        return self.sim.people.uid.isin(new) # As a boolean array

    def clear_infection(self, uids):
        """ Cure people (e.g. treated): they are susceptible again, and their serology stays as it is """
        self.move(uids, self.susceptible)
        return

    def init_post(self):
        """ Those infected at the start have the serology and symptoms of a new infection """
        super().init_post()
        self.set_infection(self.infected.uids)
        return

    def set_infection(self, uids):
        """ At infection: serology to come, and whether the chancre will be visible """
        p = self.pars
        self.ever_exposed[uids] = True
        self.n_infections[uids] += 1

        # Treponemal seroconversion after the window period; 80% stay positive for life, the others sero-revert later
        self.ti_trep_start[uids] = self.ti + rr(p.dur_to_trep.rvs(uids))
        reverters = uids[~p.p_trep_persists.rvs(uids)]
        self.ti_trep_end[reverters] = self.ti_trep_start[reverters] + rr(p.time_to_trep_revert.rvs(reverters))

        # Whether the primary chancre will be at a visible site
        p.p_symp_primary_dist.set(self.by_sex(p.p_symp_primary, uids))
        self.chancre_visible[uids] = p.p_symp_primary_dist.rvs(uids)
        return

    def step_state(self):
        """ Take the declared transitions; then update serology, visible symptoms, and the infectiousness of latent syphilis """
        super().step_state()
        ti = self.ti
        p = self.pars

        # Entering the secondary stage (from primary, or a relapse): the non-treponemal titre rises, and the rash may be noticed
        secondary = self.secondary_onset.uids | self.reactivation.uids
        self.nontrep[secondary] = True
        p.p_symp_secondary_dist.set(self.by_sex(p.p_symp_secondary, secondary))
        self.rash_visible[secondary] = p.p_symp_secondary_dist.rvs(secondary)

        # Entering the latent stage: the rash clears, and infectiousness begins to decay
        latent = self.latency.uids
        self.rash_visible[latent] = False
        self.ti_latent[latent] = ti

        # Entering late latent: the non-treponemal titre falls below detection some years later
        late = self.to_late.uids
        self.ti_nontrep_end[late] = ti + rr(p.time_to_undetectable.rvs(late))

        # Serology over time
        self.nontrep[self.nontrep & (self.ti_nontrep_end <= ti)] = False
        self.trep[~self.trep & (self.ti_trep_start <= ti)] = True
        self.trep[self.trep & (self.ti_trep_end <= ti)] = False

        # Susceptibility and infectiousness of each person (which connectors change after this), with the decay in the latent stage
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1
        latent = self.latent.uids
        half_life = p.rel_trans_latent_half_life.to_dt()
        self.rel_trans[latent] = 0.5**((ti - self.ti_latent[latent])/half_life)
        return

    def step(self):
        """ Transmission; then the people infected have the serology and symptoms of a new infection """
        super().step()
        self.set_infection(self.infection.uids)
        return

    def init_results(self):
        """ Initialize results """
        super().init_results()
        self.define_sti_results(include_care=self.include_care)
        results = [
            ss.Result('cum_infections', dtype=int, label='Cumulative infections'),
            ss.Result('trep_prevalence', dtype=float, scale=False, label="Treponemal-test positive (ever-exposed proxy)"),
            ss.Result('pregnant_prevalence', dtype=float, scale=False, label="Pregnant prevalence", auto_plot=False),
            ss.Result('detected_pregnant_prevalence', dtype=float, scale=False, label="ANC prevalence", auto_plot=False),
            ss.Result('delivery_prevalence', dtype=float, scale=False, label="Delivery prevalence", auto_plot=False),
            ss.Result('sexually_transmissible_prevalence', dtype=float, scale=False, label="Sexually transmissible (primary + secondary + early latent)"),
            ss.Result('symptomatic_prevalence', dtype=float, scale=False, label="Symptomatic stage prevalence (primary + secondary)"),
            ss.Result('primary_prevalence', dtype=float, scale=False, label="Primary stage prevalence"),
            ss.Result('nontrep_prevalence', dtype=float, scale=False, label="Non-treponemal (RPR) positive"),
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

        self.define_results(*results)
        return

    def update_results(self):
        super().update_results()
        self.update_sti_results()
        ti = self.ti
        res = self.results
        ppl = self.sim.people
        new_inf = self.new_infected()

        adults = (ppl.age >= self.age_range[0]) & (ppl.age < self.age_range[1])
        sexually_active_adults = adults & self.sim.networks.structuredsexual.active(self.sim.people)

        # Prevalence among sexually active adults
        res['prevalence'][ti] = cond_prob(self.infected, sexually_active_adults)
        res['trep_prevalence'][ti] = cond_prob(self.trep, sexually_active_adults)
        res['nontrep_prevalence'][ti] = cond_prob(self.nontrep, sexually_active_adults)
        res['sexually_transmissible_prevalence'][ti] = cond_prob(self.sexually_transmissible, sexually_active_adults)
        res['symptomatic_prevalence'][ti] = cond_prob(self.symptomatic, sexually_active_adults)
        res['primary_prevalence'][ti] = cond_prob(self.primary, sexually_active_adults)

        # Pregnant women, if there is pregnancy
        if 'pregnancy' in self.sim.demographics.keys():
            preg_prev = cond_prob(self.infected, ppl.pregnancy.P)
            res['pregnant_prevalence'][ti] = preg_prev
            res['detected_pregnant_prevalence'][ti] = preg_prev * self.pars.anc_detection
            res['delivery_prevalence'][ti] = cond_prob(self.infected, ppl.pregnancy.ti_delivered == ti)

        res['new_reinfections'][ti] = count(new_inf & (self.n_infections > 1))

        # Add FSW and clients to results:
        if self.store_sw:
            fsw = self.sim.networks.structuredsexual.fsw
            clients = self.sim.networks.structuredsexual.client
            res['prevalence_sw'][ti] = cond_prob(self.infected, fsw)
            res['new_infections_sw'][ti] = count(new_inf & fsw)
            res['new_infections_not_sw'][ti] = count(new_inf & ~fsw)
            res['prevalence_client'][ti] = cond_prob(self.infected, clients)
            res['new_infections_client'][ti] = count(new_inf & clients)
            res['new_infections_not_client'][ti] = count(new_inf & ~clients)

        # Add risk groups
        if self.store_risk_groups:
            for risk_group in range(self.sim.networks.structuredsexual.pars.n_risk_groups):
                for sex in ['female', 'male']:
                    prev_denom = (self.sim.networks.structuredsexual.risk_group == risk_group) & (self.sim.people[sex])
                    risk_group_new_inf = new_inf & (self.sim.networks.structuredsexual.risk_group == risk_group) & (self.sim.people[sex])
                    if risk_group_new_inf.any():
                        res['prevalence_risk_group_' + str(risk_group) + '_' + sex][ti] = cond_prob(self.infected, prev_denom)
                        res['new_infections_risk_group_' + str(risk_group) + '_' + sex][ti] = count(risk_group_new_inf)

        # Results by age and sex
        for skk, skl, attr in self.sex_strata():
            sex_denom = sexually_active_adults & ppl[attr]
            if skk != '':
                res[f'prevalence{skk}'][ti] = cond_prob(self.infected, sex_denom)
                res[f'trep_prevalence{skk}'][ti] = cond_prob(self.trep, sex_denom)
                res[f'nontrep_prevalence{skk}'][ti] = cond_prob(self.nontrep, sex_denom)
                res[f'symptomatic_prevalence{skk}'][ti] = cond_prob(self.symptomatic, sex_denom)
                res[f'sexually_transmissible_prevalence{skk}'][ti] = cond_prob(self.sexually_transmissible, sex_denom)
                res[f'primary_prevalence{skk}'][ti] = cond_prob(self.primary, sex_denom)

            # Age-binned prevalences (denominator = all people of that sex in each age bin)
            if self.age_bins is not None:
                ppl_hist = self.agehist(ppl[attr])
                age_results = dict(
                    sexually_transmissible_prevalence = div(self.agehist(self.sexually_transmissible & ppl[attr]), ppl_hist),
                    symptomatic_prevalence = div(self.agehist(self.symptomatic & ppl[attr]), ppl_hist),
                    primary_prevalence = div(self.agehist(self.primary & ppl[attr]), ppl_hist),
                    trep_prevalence = div(self.agehist(self.trep & ppl[attr]), ppl_hist),
                    nontrep_prevalence = div(self.agehist(self.nontrep & ppl[attr]), ppl_hist),
                )
                for akey, ares in age_results.items():
                    for ai, (ab1, ab2) in enumerate(self.age_strata()):
                        res[f'{akey}{skk}_{ab1}_{ab2}'][ti] = ares[ai]
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        res['cum_infections'][:] = np.cumsum(res['new_infections'][:])
        res['cum_congenital'][:] = np.cumsum(res['new_congenital'])
        res['cum_congenital_deaths'][:] = np.cumsum(res['new_congenital_deaths'])
        return
