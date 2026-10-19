"""
Define syphilis disease module
"""

import numpy as np
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI, sexed

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


# The congenital outcomes, and the mother's stages they depend on
outcome_keys = ['miscarriage', 'nnd', 'stillborn', 'congenital', 'normal']
outcome_stages = ['active', 'early', 'late']
v3_stages = dict(mat_active='active', early='early', late='late')


def birth_outcome_pars(birth_outcomes):
    """
    The probabilities of each congenital outcome by the mother's stage, from STIsim 3's `birth_outcomes` (a dict by
    stage, `mat_active`, `early`, `late`, of `ss.choice(a=5, p=[miscarriage, nnd, stillborn, congenital, normal])` or
    of the list of probabilities), as the parameters `p_<outcome>_<stage>` of `sti.Syphilis`
    """
    out = {}
    for key,val in birth_outcomes.items():
        stage = v3_stages.get(key, key)
        if stage not in outcome_stages:
            errormsg = f'birth_outcomes is by stage ({sc.strjoin(v3_stages.keys())}), not {key!r}'
            raise ValueError(errormsg)
        p = np.asarray(val.pars.p if isinstance(val, ss.Dist) else val, dtype=float)
        if len(p) != 5 or not np.isclose(p.sum(), 1):
            errormsg = f'birth_outcomes[{key!r}] must be 5 probabilities ({sc.strjoin(outcome_keys)}) that sum to 1, not {p}'
            raise ValueError(errormsg)
        for k,v in zip(outcome_keys[:4], p[:4]):
            out[f'p_{k}_{stage}'] = float(v)
    return out


class Syphilis(BaseSTI):
    """
    Syphilis: adult stages, serology, and congenital syphilis by the mother's stage

    States: S (never infected), Sp and Spl (cured in an early or a late stage: treponemal
    test positive), Sr (sero-reverted: treponemal negative); E (incubating), Pv and Ph
    (primary, with the chancre at a visible or hidden site), S2v and S2h (secondary, rash
    visible or not), EL and LL (early and late latent), T (tertiary); and, for fetuses,
    C<outcome> (infected in utero, by outcome) and CS (born with congenital syphilis).

        S, Sp, Spl, Sr -> E                     sexual transmission
        E -> Pv (p_symp_primary) | Ph           after dur_exposed
        Pv, Ph -> S2v (p_symp_secondary) | S2h  after dur_primary; RPR positive
        S2v, S2h -> EL                          after dur_secondary
        EL -> S2v | S2h                         p_reactivate, after time_to_reactivate (reactivate_v, _h); RPR positive
        EL -> LL                                the rest, after dur_early
        LL -> T                                 p_tertiary, after time_to_tertiary from latent onset; the rest stay
        T -> dead                               p_death, after time_to_death; the rest stay
        E, Pv, Ph, S2v, S2h, EL -> Sp           cure_early: moved only by treatment
        LL, T -> Spl                            cure_late: moved only by treatment
        Sp, Spl -> Sr                           1 - p_trep_persists, after time_to_trep_revert; the rest stay

    Infectiousness is `rel_trans_primary` in primary, `rel_trans_secondary` in secondary
    and at the start of latency, halving every `rel_trans_latent_half_life` of time since
    latent onset (`(EL | LL).time`), and `rel_trans_tertiary` in tertiary; E does not
    transmit. The non-treponemal (RPR) titre is the `rpr` dimension: positive from
    secondary, falling below detection after a mean `time_to_undetectable` in late latent
    (and after late treatment), and after a mean `dur_nontrep_revert` after early treatment.

    A fetus of a mother in P, S2, EL, or LL is infected at `beta_m2c` per month, and its
    outcome is drawn from her stage (active: primary or secondary; early or late latent):
    a miscarriage or a stillbirth at delivery, a neonatal death or congenital syphilis at
    birth, or a normal birth (`p_<outcome>_<stage>`, or STIsim 3's `birth_outcomes`).
    Treating the mother can cure the fetus (`fetal_cure`), and a newborn with congenital
    syphilis can be treated (`cs_cure`).

    Results: the counts and flows of every state and transition; `new_infections`,
    `new_reinfections`, `n_infected`, `prevalence` (among adults in `age_range`),
    `n_active`, `symptomatic_prevalence`, `primary_prevalence`,
    `sexually_transmissible_prevalence`, `trep_prevalence`, `nontrep_prevalence`,
    `new_nnds`, `new_stillborns`, `new_miscarriages`, `new_congenital`,
    `new_congenital_deaths`, `new_deaths`; by sex (`_f`, `_m`) and with `.select()`/`.by()`
    by sex, age, and `rpr`.

    Args:
        init_prev (float): initial prevalence of new infections (E)
        init_latent_prev (float): initial prevalence of latent infection, spread over latency
        birth_outcomes (dict): STIsim 3's congenital outcomes by the mother's stage (see `birth_outcome_pars()`)
        kwargs (dict): parameter values by name, and those of `sti.BaseSTI`

    Examples:
        ```python
        syph = sti.Syphilis(beta_m2f=0.2, beta_m2c=0.5, init_prev=0.02, init_latent_prev=0.02)
        sim = sti.Sim(diseases=syph, demographics='zimbabwe')
        ```
    """
    S, Sp, Spl, Sr = ss.states('S Sp Spl Sr')
    E, Pv, Ph, S2v, S2h, EL, LL, T = ss.states('E Pv Ph S2v S2h EL LL T')
    Cmisc, Cnnd, Cstill, Ccs, Cnorm, CS = ss.states('Cmisc Cnnd Cstill Ccs Cnorm CS')

    # Adult natural history
    dur_exposed = ss.par(ss.normal(ss.days(50), ss.days(10)), doc='Incubation period (https://pubmed.ncbi.nlm.nih.gov/9101629/)')
    p_symp_primary = ss.par(sexed([0.3, 0.5]), doc='Probability that the primary chancre is at a visible site (women, men)')
    dur_primary = ss.par(ss.normal(ss.weeks(6), ss.weeks(1)), doc='Duration of primary syphilis (https://pubmed.ncbi.nlm.nih.gov/9101629/)')
    p_symp_secondary = ss.par(sexed([0.85, 0.85]), doc='Probability that the secondary rash is noticeable (women, men)')
    dur_secondary = ss.par(ss.lognorm_ex(ss.months(3.6), ss.months(1.5)), doc='Duration of secondary syphilis (https://pubmed.ncbi.nlm.nih.gov/9101629/)')
    dur_early = ss.par(ss.uniform(ss.months(22), ss.months(24)), doc='Duration of early latent syphilis, from latent onset')
    p_reactivate = ss.par(0.35, doc='Probability of reactivation from latent to secondary')
    time_to_reactivate = ss.par(ss.lognorm_ex(ss.years(1), ss.years(1)), doc='Time from latent onset to reactivation')
    p_tertiary = ss.par(0.35, doc='Probability of tertiary syphilis (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917057/)')
    time_to_tertiary = ss.par(ss.normal(ss.years(20) - ss.months(23), ss.years(2)), doc='Time from late latent onset to tertiary (20 years from latent onset)')
    p_death = ss.par(0.05, doc='Probability of death from tertiary syphilis')
    time_to_death = ss.par(ss.lognorm_ex(ss.years(5), ss.years(5)), doc='Time from tertiary to death')

    # Serology
    time_to_undetectable = ss.par(ss.years(15), doc='Mean time in late latent until the non-treponemal (RPR) titre falls below detection')
    dur_nontrep_revert = ss.par(ss.months(9), doc='Mean time after early treatment until the RPR titre falls below detection')
    p_trep_persists = ss.par(0.80, doc='Probability that a treponemal test stays positive for life after cure')
    time_to_trep_revert = ss.par(ss.normal(ss.years(10), ss.years(5)), doc='Time to treponemal sero-reversion, for those who revert')

    # Transmission
    rel_trans_primary = ss.par(1.0, doc='Relative infectiousness in primary syphilis')
    rel_trans_secondary = ss.par(1.0, doc='Relative infectiousness in secondary syphilis, and at latent onset')
    rel_trans_tertiary = ss.par(0.0, doc='Relative infectiousness in tertiary syphilis')
    rel_trans_latent_half_life = ss.par(ss.years(1), doc='Half-life of infectiousness in latent syphilis')

    # Congenital outcomes by the mother's stage (normal is the rest)
    p_miscarriage_active = ss.par(0.00, doc='Probability of miscarriage, mother in primary or secondary')
    p_nnd_active         = ss.par(0.25, doc='Probability of neonatal death, mother in primary or secondary')
    p_stillborn_active   = ss.par(0.25, doc='Probability of stillbirth, mother in primary or secondary')
    p_congenital_active  = ss.par(0.45, doc='Probability of congenital syphilis, mother in primary or secondary')
    p_miscarriage_early  = ss.par(0.00, doc='Probability of miscarriage, mother in early latent')
    p_nnd_early          = ss.par(0.05, doc='Probability of neonatal death, mother in early latent')
    p_stillborn_early    = ss.par(0.10, doc='Probability of stillbirth, mother in early latent')
    p_congenital_early   = ss.par(0.40, doc='Probability of congenital syphilis, mother in early latent')
    p_miscarriage_late   = ss.par(0.00, doc='Probability of miscarriage, mother in late latent')
    p_nnd_late           = ss.par(0.00, doc='Probability of neonatal death, mother in late latent')
    p_stillborn_late     = ss.par(0.10, doc='Probability of stillbirth, mother in late latent')
    p_congenital_late    = ss.par(0.10, doc='Probability of congenital syphilis, mother in late latent')
    mtct_rate = ss.par(ss.permonth(0), doc='Rate of in-utero infection (from beta_m2c)')

    # Sexual transmission: by stage; latent infectiousness decays (see modifiers())
    infection   = S >> E @ BaseSTI.beta*(Pv + Ph + S2v + S2h + EL + LL + T)
    reinfection = (Sp | Spl | Sr) >> E @ BaseSTI.beta*(Pv + Ph + S2v + S2h + EL + LL + T)

    # Natural history
    primary    = ss.Transition(E, [Pv, Ph], dur=dur_exposed, p=[p_symp_primary, None])
    secondary  = ss.Transition([Pv, Ph], [S2v, S2h], dur=dur_primary, p=[p_symp_secondary, None], move=dict(rpr=True))
    latent     = ss.Transition([S2v, S2h], EL, dur=dur_secondary)
    reactivate_v = ss.Transition(EL, S2v, dur=time_to_reactivate, p=p_reactivate*p_symp_secondary, move=dict(rpr=True)) # Split at latent onset
    reactivate_h = ss.Transition(EL, S2h, dur=time_to_reactivate, p=p_reactivate*(1 - p_symp_secondary), move=dict(rpr=True))
    to_late    = ss.Transition(EL, LL, dur=dur_early, p=None) # The rest: those who will not reactivate
    tertiary   = ss.Transition(LL, T, dur=time_to_tertiary, p=p_tertiary) # The rest stay
    death      = ss.Transition(T, ss.null, dur=time_to_death, p=p_death) # The rest stay
    trep_revert = ss.Transition([Sp, Spl], Sr, dur=time_to_trep_revert, p=1 - p_trep_persists) # The rest stay

    # Treatment
    cure_early = ss.Transition(E | Pv | Ph | S2v | S2h | EL, Sp, undriven=True)
    cure_late  = ss.Transition(LL | T, Spl, undriven=True)

    # Congenital: the outcome is drawn from the mother's stage when the fetus is infected
    mtct = ss.Transition(S, [Cmisc, Cnnd, Cstill, Ccs, Cnorm], rate=mtct_rate*ss.mother(Pv | Ph | S2v | S2h | EL | LL), p=[
        p_miscarriage_active*ss.mother(Pv | Ph | S2v | S2h) + p_miscarriage_early*ss.mother(EL) + p_miscarriage_late*ss.mother(LL),
        p_nnd_active*ss.mother(Pv | Ph | S2v | S2h) + p_nnd_early*ss.mother(EL) + p_nnd_late*ss.mother(LL),
        p_stillborn_active*ss.mother(Pv | Ph | S2v | S2h) + p_stillborn_early*ss.mother(EL) + p_stillborn_late*ss.mother(LL),
        p_congenital_active*ss.mother(Pv | Ph | S2v | S2h) + p_congenital_early*ss.mother(EL) + p_congenital_late*ss.mother(LL),
        None])
    miscarriage = ss.Transition(Cmisc, ss.null, on=ss.state('pregnancy.delivery')) # A fetal loss: classified by the pregnancy
    stillbirth  = ss.Transition(Cstill, ss.null, on=ss.state('pregnancy.delivery'))
    nnd         = ss.Transition(Cnnd, ss.null, on=ss.state('pregnancy.birth'))
    born_cs     = ss.Transition(Ccs, CS, on=ss.state('pregnancy.birth'))
    born_ok     = ss.Transition(Cnorm, S, on=ss.state('pregnancy.birth')) # Susceptible later in life (in STIsim 3 they never were)
    fetal_cure  = ss.Transition(Cmisc | Cnnd | Cstill | Ccs, Cnorm, undriven=True) # Treatment of the mother
    cs_cure     = ss.Transition(CS, S, undriven=True) # Treatment of the newborn

    # Results
    n_infected = E + Pv + Ph + S2v + S2h + EL + LL + T
    prevalence = (E + Pv + Ph + S2v + S2h + EL + LL + T)/(S + Sp + Spl + Sr + E + Pv + Ph + S2v + S2h + EL + LL + T)
    incidence = infection/S
    new_infections = infection + reinfection
    new_reinfections = ss.outcome(reinfection)
    n_active = Pv + Ph + S2v + S2h
    symptomatic_prevalence = (Pv + Ph + S2v + S2h)/(S + Sp + Spl + Sr + E + Pv + Ph + S2v + S2h + EL + LL + T)
    primary_prevalence = (Pv + Ph)/(S + Sp + Spl + Sr + E + Pv + Ph + S2v + S2h + EL + LL + T)
    sexually_transmissible_prevalence = (Pv + Ph + S2v + S2h + EL)/(S + Sp + Spl + Sr + E + Pv + Ph + S2v + S2h + EL + LL + T)
    trep_prevalence = (Pv + Ph + S2v + S2h + EL + LL + T + Sp + Spl)/(S + Sp + Spl + Sr + E + Pv + Ph + S2v + S2h + EL + LL + T)
    new_nnds = ss.outcome(nnd)
    new_stillborns = ss.outcome(stillbirth)
    new_miscarriages = ss.outcome(miscarriage)
    new_congenital = ss.outcome(born_cs)
    new_congenital_deaths = nnd + stillbirth + miscarriage
    new_deaths = ss.outcome(death)
    reactivations = reactivate_v + reactivate_h

    defaults = dict(beta_m2f=0.1, eff_condom=0.0, init_prev=0)

    def __init__(self, name='syph', init_latent_prev=None, birth_outcomes=None, init_prev_data=None, init_prev_latent_data=None, **kwargs):
        if init_prev_latent_data is not None:
            errormsg = 'init_prev_latent_data (initial latent prevalence by risk group) is not supported yet in Starsim 4; give init_latent_prev'
            raise NotImplementedError(errormsg)
        for key in ['dist_ti_init_infected', 'rel_init_prev', 'anc_detection', 'dur_to_trep']: # Starsim 3 parameters with no Starsim 4 equivalent
            kwargs.pop(key, None)
        pars = kwargs.pop('pars', None) or {}
        init_latent_prev = sexed(pars.pop('init_latent_prev', init_latent_prev))
        birth_outcomes = pars.pop('birth_outcomes', birth_outcomes)
        if birth_outcomes is not None:
            kwargs.update(birth_outcome_pars(birth_outcomes))
        for key in ['p_reactivate', 'p_tertiary', 'p_death', 'p_trep_persists']: # Starsim 3: bernoullis
            for src in [pars, kwargs]:
                if isinstance(src.get(key), ss.bernoulli):
                    src[key] = src[key].pars.p
        if init_latent_prev:
            kwargs['init'] = sc.mergedicts(kwargs.get('init'), dict(EL=ss.spread(init_latent_prev)))
        super().__init__(name=name, pars=pars, init_prev_data=init_prev_data, **kwargs)
        st = self.graph.states
        self.latent_time = (st.EL | st.LL).time # Made before the sim, so the time of entering latency is kept
        self.stratify(rpr=[False, True], split=[1, 0],
                      rpr_late=ss.transfer(rpr=False, dur=self.pars.time_to_undetectable, states=['LL', 'T', 'Spl', 'Sr']),
                      rpr_tx=ss.transfer(rpr=False, dur=self.pars.dur_nontrep_revert, states=['Sp']))
        return

    def betamap(self, sim):
        """ As for other STIs, but nothing over the maternal network: in-utero infection is `mtct`, from the mother's stage """
        out = super().betamap(sim)
        for name in out:
            if ss.standardize_netkey(name) == 'maternal':
                out[name] = 0
        return out

    def init_pre(self, sim):
        m2c = self.value('beta_m2c')
        self.pars.mtct_rate.value = ss.permonth(m2c or 0)
        super().init_pre(sim)
        return

    def modifiers(self, sim):
        """ Infectiousness by stage; in latency it halves every rel_trans_latent_half_life since latent onset """
        st, p, name = self.graph.states, self.pars, self.name
        latent = st.EL | st.LL
        rt = f'{name}.rel_trans'
        return [
            ss.modify(rt, multiply=p.rel_trans_primary, where=st.Pv | st.Ph, name=f'{name}_rel_trans_primary'),
            ss.modify(rt, multiply=p.rel_trans_secondary, where=st.S2v | st.S2h | latent, name=f'{name}_rel_trans_secondary'),
            ss.modify(rt, multiply=p.rel_trans_tertiary, where=st.T, name=f'{name}_rel_trans_tertiary'),
            ss.modify(rt, multiply=0.5**(self.latent_time/p.rel_trans_latent_half_life), where=latent, name=f'{name}_latent_decay'),
        ]

    def flat_results(self):
        return super().flat_results() + [
            ('prevalence', 'prevalence', True, 'Prevalence'),
            ('n_active', 'n_active', None, 'Number with primary or secondary syphilis'),
            ('symptomatic_prevalence', 'symptomatic_prevalence', True, 'Symptomatic prevalence'),
            ('primary_prevalence', 'primary_prevalence', True, 'Primary stage prevalence'),
            ('sexually_transmissible_prevalence', 'sexually_transmissible_prevalence', True, 'Sexually transmissible prevalence'),
            ('trep_prevalence', 'trep_prevalence', True, 'Treponemal-test positive'),
        ]

    def finalize_results(self):
        """ Add the non-treponemal (RPR) prevalence and the cumulative congenital results """
        super().finalize_results()
        res = self.results
        states = ['S', 'Sp', 'Spl', 'Sr', 'E', 'Pv', 'Ph', 'S2v', 'S2h', 'EL', 'LL', 'T']
        pos = sum(res[f'n_{s}'].by('rpr').values[:, 1] for s in states)
        tot = sum(res[f'n_{s}'].values for s in states)
        vals = np.divide(pos, tot, out=np.full(len(tot), np.nan), where=tot > 0)
        tv = res.n_S.timevec
        self.results += ss.Result('nontrep_prevalence', values=vals, scale=False, label='Non-treponemal (RPR) positive', timevec=tv, module=self.name)
        for key in ['congenital', 'congenital_deaths']:
            self.results += ss.Result(f'cum_{key}', values=np.cumsum(res[f'new_{key}'].values), scale=False, label=f'Cumulative {key.replace("_", " ")}', timevec=tv, auto_plot=False, module=self.name)
        return

    # States as STIsim 3 named them, as agents, for connectors and custom code
    @property
    def susceptible(self): return self.S | self.Sp | self.Spl | self.Sr
    @property
    def exposed(self):     return self.E
    @property
    def primary_(self):    return self.Pv | self.Ph
    @property
    def secondary_(self):  return self.S2v | self.S2h
    @property
    def early(self):       return self.EL
    @property
    def late(self):        return self.LL
    @property
    def latent_(self):     return self.EL | self.LL
    @property
    def tertiary_(self):   return self.T
    @property
    def symptomatic(self): return self.Pv | self.Ph | self.S2v | self.S2h
    @property
    def infectious(self):  return self.Pv | self.Ph | self.S2v | self.S2h | self.EL | self.LL
    @property
    def infected(self):    return self.E | self.infectious | self.T
    @property
    def mat_active(self):  return self.E | self.Pv | self.Ph | self.S2v | self.S2h
    @property
    def sexually_transmissible(self): return self.Pv | self.Ph | self.S2v | self.S2h | self.EL
    @property
    def ulcerative(self):  return self.Pv
    @property
    def visibly_symptomatic(self): return self.Pv | self.S2v
    @property
    def sus_not_naive(self): return self.Sp | self.Spl | self.Sr
    @property
    def congenital(self):  return self.CS
    @property
    def treatable(self):   return self.infected
