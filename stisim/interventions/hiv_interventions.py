"""
HIV testing, treatment, and prevention

- `HIVTest` (and `InfantHIVTest`): who is diagnosed, as HIV's `diagnosed`, with a referral to ART
- `ART`: antiretroviral therapy, a declared module of its own (naive, on ART, off ART), whose effects on HIV are modifiers
- `VMMC` and `Prep`: products (`ss.Vaccine`) delivered to a coverage target
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
from stisim.interventions.base_interventions import STITest
from stisim.interventions.utils import parse_coverage, age_sex_mask

__all__ = ["HIVDx", "HIVTest", "InfantHIVTest", "ART", "VMMC", "Prep"]


def raw(arr):
    """ A boolean array over all agents (alive or not), from a BoolArr or an array """
    return arr.raw if hasattr(arr, 'raw') else np.asarray(arr)


def find_mod(sim, name):
    """ A module of the sim by name, among its interventions and diseases """
    return sim.interventions[name] if name in sim.interventions else sim.diseases[name]


def coverage_keys(spec, valid_names, sex=None):
    """ The strata of a coverage target, without a sim: [] for none, [None] for one target, or a list of (age bin, sex) """
    cov, fmt, age_bins, sex_keys = parse_coverage(spec, valid_names=valid_names, yearvec=np.array([2000.0]))
    if cov is None:
        return []
    if not isinstance(cov, dict):
        return [None]
    return [(ab, sk if sk is not None else sex) for ab in age_bins for sk in (sex_keys or [None])]


def find_art(sim):
    """ The ART module of a sim, or None """
    for mod in sim.interventions.values():
        if isinstance(mod, ART):
            return mod
    return None


class set_attr(ss.modify):
    """
    Replace an attribute's value while a condition holds: `ss.modify(target, set=...)` for an attribute (e.g. `hiv.cd4`)

    UPSTREAM: ss.modify() refuses set= on anything but a parameter (starsim/connectors.py:332), though the
    agent engine applies set= to an attribute (Agents.modify_attr()); this makes the modifier as a multiplier
    and switches it to set= once it is checked.
    """
    def __init__(self, target, set, **kwargs):
        super().__init__(target, multiply=set, **kwargs)
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.op = 'set'
        return


def add_modules(sim, connectors=(), interventions=()):
    """
    Add modules that a module brings with it (modifiers and deliveries) to the sim, during its init_pre()

    UPSTREAM: a module cannot declare the modifiers and deliveries that come with it (BaseSTI.add_modifiers has the same TODO)
    """
    for kind, mods in [('connectors', connectors), ('interventions', interventions)]:
        for mod in mods:
            if mod.name not in sim[kind]:
                sim[kind][mod.name] = mod
                mod.init_pre(sim)
    return


def coverage_strata(spec, valid_names, sim, missing_fill=0.0, smoothness=0, format_priority='n', sex=None):
    """
    A coverage target in any format `parse_coverage()` accepts, as a list of (stratum, share) pairs: the stratum is
    None or a function of the sim giving who is in it, and the share a function of the sim and of the number of
    people it is a share of, giving the share to cover now (converting a number of people, 'n', to a share)

    Args:
        sex (int): for age-only strata, the sex they are of (e.g. 1 for VMMC); None for both
    """
    cov, fmt, age_bins, sex_keys = parse_coverage(spec, valid_names=valid_names, yearvec=sim.t.yearvec,
                                                  smoothness=smoothness, format_priority=format_priority, missing_fill=missing_fill)
    if cov is None:
        return [], (cov, fmt, age_bins, sex_keys)

    def share_of(arr):
        def share(sim, n_people):
            ti = sim.ti
            val = arr[ti] if ti < len(arr) else arr[-1]
            f = fmt[ti] if isinstance(fmt, np.ndarray) else fmt
            if f == 'n':
                return float(np.clip(sc.safedivide(val, n_people), 0, 1))
            return float(np.clip(val, 0, 1))
        return share

    if not isinstance(cov, dict):
        return [(None, share_of(cov))], (cov, fmt, age_bins, sex_keys)
    out = []
    for ab in age_bins:
        for sk in (sex_keys or [None]):
            key = (ab, sk) if sk is not None else ab
            arr = cov.get(key, np.full(len(sim.t.yearvec), missing_fill))
            s = sk if sk is not None else sex
            stratum = (lambda ab, s: lambda sim: age_sex_mask(ab, s, sim.people))(ab, s)
            stratum.label = f'{ab}' + ('' if s is None else f' {["f", "m"][s]}')
            out.append((stratum, share_of(arr)))
    return out, (cov, fmt, age_bins, sex_keys)


def stock_share(share, among):
    """ An `ss.share` whose value is worked out from the number of people it is a share of, so a number of people ('n') becomes a share """
    def value(sim):
        sel = among(sim)
        sel = sel.raw if hasattr(sel, 'raw') else np.asarray(sel)
        n_people = np.count_nonzero(sel & sim.people.alive.raw)*sim._current_scale # UPSTREAM: people per agent is private (sim._current_scale)
        return share(sim, n_people)
    return ss.share(value, among=among)


# %% Testing

class HIVDx(ss.Product):
    """
    A perfect HIV test, used by `HIVTest`: positive if infected, negative if not

    Args:
        *args: passed to `ss.Product`
        **kwargs: passed to `ss.Product`
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.result_list = ['positive', 'negative']

    def administer(self, sim, uids):
        hiv = sim.diseases.hiv
        outcomes = {r: ss.uids() for r in self.result_list}
        outcomes['positive'] = hiv.infected.uids.intersect(uids)
        outcomes['negative'] = hiv.susceptible.uids.intersect(uids)
        return outcomes


class HIVTest(STITest):
    """
    HIV testing

    Tests eligible people (by default, those not yet diagnosed). A positive result sets `hiv.diagnosed` and
    `hiv.ti_diagnosed`, and refers the person to ART (if `sti.ART` is in the sim), who is ready to start it after
    `dur_dx2tx`.

    Args:
        product: the test (default: a perfect `HIVDx`)
        test_prob_data: the probability of testing: annual if `dt_scale=True` (the default), per step otherwise; an array gives one per year in `years`
        years (array): the years of `test_prob_data` (and the start and stop)
        start (float): the year testing starts
        eligibility (func): who can be tested, a function of the sim (default: the undiagnosed)
        newborn_test (`InfantHIVTest`): if given, the children of mothers who test positive in pregnancy are tested at birth
        dur_dx2tx (ss.Dist): the time from diagnosis to being ready to start ART (default: none)
        dt_scale (bool): whether `test_prob_data` is annual (True) or per step (False)

    Examples:
        ```python
        test = sti.HIVTest(test_prob_data=0.2, start=2000, name='hiv_test')
        fsw_test = sti.HIVTest(test_prob_data=0.5, name='fsw_test',
                               eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed)
        ```
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None,
                 eligibility=None, name=None, label=None, newborn_test=None, **kwargs):
        if product is None:
            suffix = f'_{name}' if name else ''
            product = HIVDx(name=f'HIVDx{suffix}')
        super().__init__(product=product, test_prob_data=test_prob_data, years=years,
                         start=start, eligibility=eligibility, name=name, label=label)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed
        self.newborn_test = newborn_test
        self.define_pars(dur_dx2tx=ss.constant(0))
        self.update_pars(pars, **kwargs)
        return

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        diagnose(sim, pos_uids, self.ti, self.pars.dur_dx2tx)

        # Test at birth the children of mothers diagnosed in pregnancy
        if self.newborn_test is not None and len(pos_uids):
            mn = sim.get_module(ss.MaternalNet, die=False)
            pregnancy = sim.get_module(ss.Pregnancy, die=False)
            if mn is None or pregnancy is None:
                raise RuntimeError('HIVTest: cannot schedule newborn_test without ss.MaternalNet and ss.Pregnancy in the sim')
            pos_mask = np.isin(mn.p1, pos_uids)
            unborn = mn.p2[pos_mask]
            ti_births = pregnancy.ti_delivery[mn.p1[pos_mask]]
            valid = ~np.isnan(ti_births)
            if valid.any():
                self.newborn_test.schedule(unborn[valid], ti_births[valid].astype(int))
        return outcomes


def diagnose(sim, uids, ti, delay=None):
    """
    Diagnose people with HIV: set `hiv.diagnosed`, and refer them to ART (if it is in the sim), ready to start after `delay`

    Args:
        sim (Sim): the sim
        uids (uids): the people diagnosed
        ti (int): the time step of diagnosis
        delay (ss.Dist): the time from diagnosis to being ready to start ART
    """
    if not len(uids):
        return
    hiv = sim.diseases.hiv
    hiv.diagnosed[uids] = True
    hiv.ti_diagnosed[uids] = ti
    art = find_art(sim)
    if art is not None:
        steps = delay.rvs(uids, round=True) if delay is not None else 0
        art.refer(uids, ti + steps)
    return


class InfantHIVTest(HIVTest):
    """
    HIV test of infants born to mothers diagnosed in pregnancy

    Scheduled by `HIVTest` (with `newborn_test=`) or `ANCTest` when the mother tests positive; tests only at the scheduled step.

    Args:
        test_prob (float): the probability that an infant scheduled for a test is tested
    """
    def __init__(self, test_prob=1.0, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label, **kwargs)
        self.test_prob_dist = ss.bernoulli(p=test_prob)

    def check_eligibility(self):
        scheduled_now = (self.ti_scheduled == self.ti).uids
        return self.test_prob_dist.filter(scheduled_now)


# %% ART

def interval_mult(bands, values):
    """ A multiplier by value from (lo, hi, mult) bands, multiplied where they overlap, for an array of values """
    out = np.ones(len(values))
    for lo, hi, mult in bands:
        out[(values >= lo) & (values < hi)] *= mult
    return out


class ART(ss.Disease, ss.cascade):
    """
    Antiretroviral therapy: a module of its own, which changes HIV through modifiers

    States (everyone starts in `naive`):

        naive -> on_art       start: by a delivery of ART, to people diagnosed and ready (see sti.HIVTest)
        off_art -> on_art     start: restarting
        on_art -> off_art     stop, after dur_on_art (x rel_dur_on_art, x dur_on_art_trend)
        off_art -> dead       death, after the post-ART survival time

    Attributes, for each person: `cd4_preart` (CD4 when starting), `cd4_potential` (the CD4 that ART takes
    them towards, from their first start), `cd4_stop` (CD4 when stopping), `suppressed` (whether ART
    suppresses their virus, drawn at each start with p_effective_art), and `care_seeking`.

    Effects on HIV (`ss.modify`):

    - infectiousness on ART falls linearly over time_to_art_efficacy to 1 - effective_art_efficacy, or 1 - nonsupp_art_efficacy if not suppressed
    - from the first start, ART replaces HIV's natural history: HIV's progression and AIDS death are stopped
    - CD4 on ART rises from cd4_preart towards cd4_potential (logistic, at art_cd4_growth per month); after
      stopping, it falls linearly from cd4_stop to 0 over the post-ART survival time, which is
      log(cd4_preart) x cd4_stop / cd4_potential years
    - HIV's CD4-banded death rate on ART is x rel_art_mortality_effective (suppressed) or x rel_art_mortality_unsupp_m/f
      (not), x art_death_age (by age), x art_death_dur (by time on ART)
    - transmission to the unborn and breastfed children of mothers on ART is x (1 - pmtct_efficacy)

    Who starts ART:

    - with no coverage: everyone ready (diagnosed, after the referral delay), once, with probability art_initiation
    - with a coverage target (number on ART, or share of people with HIV; constant, by year, or by age and sex):
      as many as the target needs among those ready (naive or off ART), lowest CD4 over care seeking first, and
      people above the target stop ART, highest first. One delivery per stratum keeps each stratum's own target.

    `vls_coverage` (a share, by year, or by age and sex) keeps the share of people on ART who are suppressed at a
    target, by an adherence rank (corrected at the start of each step, agents only).

    Results: `n_naive`, `n_on_art`, `n_off_art`, the flows `start`, `stop`, and `death`; `p_on_art` (share of
    people with HIV), `new_agents_on_art`, `n_on_art_pregnant`, `n_vls` and `p_vls` (of people with HIV).
    HIV's `new_deaths` includes the deaths after stopping ART.

    Args:
        coverage: a coverage target (see `sti.parse_coverage`: a number, a dict, or a dataframe with `n_art`/`p_art`, or stratified)
        vls_coverage: a target for the share suppressed among people on ART (a share, a dict, or a dataframe with `p_vls`, or stratified)
        smoothness (float): interpolation of coverage over time (0 for linear)
        format_priority (str): with both `n_art` and `p_art`, which to use ('n' or 'p')
        disease (str): the name of the HIV module (default 'hiv')
        kwargs (dict): parameter values by name

    Examples:
        ```python
        art = sti.ART(coverage=0.8)
        art = sti.ART(coverage={'year': [2000, 2010, 2025], 'value': [0, 0.5, 0.9]}, dur_on_art=ss.lognorm_ex(ss.years(5), ss.years(2)))
        art = sti.ART(coverage=pd.read_csv('art_coverage.csv').set_index('year'))
        ```
    """
    has_product = False # A module of states, passed with the interventions
    naive, on_art, off_art = ss.states('naive on_art off_art')

    # Who starts, and for how long
    art_initiation = ss.par(0.9, doc='Probability that someone diagnosed and ready starts ART, with no coverage target')
    dur_on_art = ss.par(ss.lognorm_ex(ss.years(3), ss.years(1.5)), doc='Time on ART before stopping')
    rel_dur_on_art = ss.par(1.0, doc='Multiplier on the time on ART')
    p_effective_art = ss.par(1.0, doc='Probability that ART suppresses the virus, drawn at each start')

    # Effects on infectiousness
    effective_art_efficacy = ss.par(0.99, doc='Reduction in infectiousness on suppressive ART (PARTNER, PARTNER2, Opposites Attract: U=U; 0.99 leaves a residual for viral blips)')
    nonsupp_art_efficacy = ss.par(0.35, doc='Reduction in infectiousness on non-suppressive ART (Quinn, Rakai 2000, at an assumed non-suppressed viral load of log10 3.75 +/- 0.8)')
    time_to_art_efficacy = ss.par(ss.months(6), doc='Time to full efficacy against transmission (linear)')
    pmtct_efficacy = ss.par(0.96, doc='Reduction in transmission to the unborn and breastfed children of mothers on ART')

    # Effects on CD4
    art_cd4_growth = ss.par(0.1, doc='Rate of CD4 recovery on ART, per month (logistic)')
    cd4_max = ss.par(1000, doc='Highest CD4 count that ART can restore')
    cd4_healthy = ss.par(500, doc='Healthy CD4 count, which scales the gain on ART')

    # Effects on mortality: a share of HIV's off-ART death rate at the same CD4, so ART is never more dangerous than none
    rel_art_mortality_effective = ss.par(0.25, doc='Share of the CD4-banded death rate kept on suppressive ART')
    rel_art_mortality_unsupp_m = ss.par(0.7, doc='Share of the CD4-banded death rate kept on non-suppressive ART, men')
    rel_art_mortality_unsupp_f = ss.par(0.35, doc='Share of the CD4-banded death rate kept on non-suppressive ART, women')

    # Transitions
    start = (naive + off_art) >> on_art
    stop = ss.Transition(on_art, off_art, dur=dur_on_art)

    # Attributes
    cd4_preart = ss.attr(1.0, on=start) # CD4 just before starting, from hiv.cd4: set by on_start()
    cd4_potential = ss.attr(1.0, on=start) # The CD4 ART takes people towards, from their first start: set by on_start()
    record_start = ss.func(lambda sim, uids: find_art(sim).on_start(uids), on=start)
    cd4_on = cd4_preart + (cd4_potential - cd4_preart)*(2/(1 + np.exp(-art_cd4_growth*on_art.time/ss.months(1))) - 1) # Logistic recovery on ART
    cd4_stop = ss.attr(cd4_on, on=stop) # CD4 when stopping
    suppressed = ss.attr(ss.bernoulli(p=p_effective_art), on=start)
    care_seeking = ss.attr(ss.normal(1, 0.5))
    death = ss.Transition(off_art, ss.null, dur=np.maximum(np.log(cd4_preart)*cd4_stop/cd4_potential, 0.1)*ss.years(1), stages=100) # Normal-ish, with a CV of 0.1, as STIsim 3's dur_post_art

    def __init__(self, pars=None, coverage=None, vls_coverage=None, smoothness=0, format_priority='n', disease='hiv',
                 art_death_age=None, art_death_dur=None, dur_on_art_trend=None, maternal_care_scale=2, name='art', **kwargs):
        pars = sc.mergedicts(pars)
        for key in ['art_death_age', 'art_death_dur', 'dur_on_art_trend', 'maternal_care_scale']: # Options in pars, as STIsim 3 took them
            if key in pars:
                locals_ = dict(art_death_age=art_death_age, art_death_dur=art_death_dur, dur_on_art_trend=dur_on_art_trend, maternal_care_scale=maternal_care_scale)
                locals_[key] = pars.pop(key)
                art_death_age, art_death_dur, dur_on_art_trend, maternal_care_scale = locals_.values()
        if 'art_cd4_pars' in kwargs or 'art_cd4_pars' in pars:
            cd4p = kwargs.pop('art_cd4_pars', pars.pop('art_cd4_pars', {}))
            pars.update({k:v for k,v in cd4p.items() if k in ['cd4_max', 'cd4_healthy']})
        for key in ['dur_post_art', 'dur_post_art_scale_factor']:
            pars.pop(key, None); kwargs.pop(key, None) # The post-ART survival time is declared (see death)
        for src in [pars, kwargs]:
            for key in ['art_initiation', 'p_effective_art']:
                if isinstance(src.get(key), ss.bernoulli):
                    src[key] = src[key].pars.p
        if 'init_prev' in kwargs:
            raise ValueError('sti.ART has no init_prev: everyone starts naive')
        names = ss.transitions.class_members(type(self)).keys()
        pars.update({k:kwargs.pop(k) for k in list(kwargs) if k in names})
        super().__init__(name=name, pars=pars, **kwargs)
        self.disease = disease
        self.raw_coverage = coverage
        self.raw_vls_coverage = vls_coverage
        self.smoothness = smoothness
        self.format_priority = format_priority
        self.art_death_age = [(0, 25, 1.0), (25, 35, 1.10), (35, 45, 1.21), (45, 125, 1.32)] if art_death_age is None else art_death_age
        self.art_death_dur = art_death_dur # (lo, hi, mult) by days on ART
        self.dur_on_art_trend = dur_on_art_trend # E.g. sc.objdict(years=[2004, 2015, 2030], vals=[0.5, 1.0, 1.5])
        self.maternal_care_scale = maternal_care_scale
        self.define_states(
            ss.FloatArr('ti_ready', label='Time ready to start ART (after diagnosis and referral)'),
            ss.FloatArr('adherence', default=ss.random(), label='Adherence rank, for vls_coverage'),
            ss.FloatArr('potential', label='CD4 potential from the first start'),
            ss.FloatArr('initiation', default=ss.random(), label='Draw for taking up ART when first ready (art_initiation), under a coverage target'),
        )
        self.steps = self.make_steps() # As a cascade's: the sim adds them to its interventions (cascade.__init__ is not called: ART's results are its own)
        return

    # Who is on ART, as STIsim 3 named it
    @property
    def on_effective_art(self):
        return self.on_art & (self.suppressed > 0.5)

    @property
    def on_nonsuppressive_art(self):
        return self.on_art & ~(self.suppressed > 0.5)

    @property
    def art_naive(self):
        return self.naive

    @property
    def art_discontinued(self):
        return self.off_art

    @property
    def ti_art(self):
        return self.ti_on_art

    def refer(self, uids, ti):
        """ Refer people who have been diagnosed to ART: ready to start it from time step `ti` (unless on ART already) """
        uids = ss.uids(uids)
        if np.ndim(ti):
            ti = np.asarray(ti)[~self.on_art[uids]]
        uids = uids[~self.on_art[uids]]
        self.ti_ready[uids] = ti
        return

    def cd4_natural(self, uids):
        """ HIV's CD4 count of some people without ART's modifiers: for people starting ART for the first time """
        hiv = self.sim.diseases[self.disease]
        key = (hiv.name, 'attr', 'cd4')
        mods = self.sim.modifiers.get(key, [])
        self.sim.modifiers[key] = [m for m in mods if not m.name.startswith(f'{self.name}_')] # UPSTREAM: see cd4_now()
        try:
            return hiv.cd4_values(uids)
        finally:
            self.sim.modifiers[key] = mods

    def cd4_now(self, uids, first=None):
        """
        HIV's CD4 count of people starting or stopping ART: for a first start, without ART; otherwise as at the end of the last step (with ART's modifiers)

        UPSTREAM: ss.attr(hiv.cd4, on=art.start) would record it at the move itself, but another module's attribute read during a move
        sees the mover in its new state (Agents._leaving covers only the module's own states), so ART's own modifiers are misapplied;
        and an attribute of another module in the class body loses its module when the module is made.
        """
        hiv = self.sim.diseases[self.disease]
        vals = np.array(hiv.cd4.raw[uids], dtype=float) # Starsim keeps derived attributes as at the end of each step
        redo = np.isnan(vals) if first is None else (first | np.isnan(vals))
        if redo.any():
            vals[redo] = self.cd4_natural(uids[redo])
        return np.clip(vals, 1, self.value('cd4_max') - 1)

    def on_start(self, uids):
        """ When people start ART: their CD4 before starting, and (the first time) the CD4 that ART takes them towards """
        is_first = np.isnan(self.potential.raw[uids])
        self.cd4_preart.raw[uids] = self.cd4_now(uids, first=is_first)
        first = uids[is_first]
        pre = self.cd4_preart.raw[first]
        cd4_max, cd4_healthy = self.value('cd4_max'), self.value('cd4_healthy')
        self.potential.raw[first] = pre*(1 + (cd4_max - pre)/cd4_healthy*np.log(cd4_max/pre))
        self.cd4_potential.raw[uids] = self.potential.raw[uids]
        return

    def check_mortality(self):
        """ On ART, the death rate must never exceed the off-ART rate at the same CD4: the largest share kept, times the largest multipliers, must be at most 1 """
        v = self.value
        age_mult = max(m for *_, m in self.art_death_age) if self.art_death_age else 1.0
        dur_mult = max(m for *_, m in self.art_death_dur) if self.art_death_dur else 1.0
        worst = max(v('rel_art_mortality_effective'), v('rel_art_mortality_unsupp_m'), v('rel_art_mortality_unsupp_f'))*age_mult*dur_mult
        if worst > 1:
            errormsg = (f'On-ART mortality can exceed off-ART mortality at the same CD4 count: max(rel_art_mortality_effective, '
                        f'rel_art_mortality_unsupp_m, rel_art_mortality_unsupp_f) * max(art_death_age mult) * max(art_death_dur mult) = {worst:.3f} > 1. '
                        'Reduce these pars so the product stays <= 1.')
            raise ValueError(errormsg)
        return

    def init_pre(self, sim):
        self.check_mortality()
        intvs = list(sim.pars.interventions.values()) if isinstance(sim.pars.interventions, dict) else sc.tolist(sim.pars.interventions)
        if not any(isinstance(m, HIVTest) for m in intvs):
            ss.warn('ART added without an HIV testing intervention; diagnosed agents may never be identified')
        hiv = sim.diseases[self.disease]
        super().init_pre(sim)
        strata, parsed = coverage_strata(self.raw_coverage, ['n_art', 'p_art'], sim, smoothness=self.smoothness, format_priority=self.format_priority)
        self.shares = [share for stratum,share in strata]
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parsed
        add_modules(sim, connectors=self.make_modifiers(sim, hiv))
        hiv.other_deaths.append((self.name, 'death')) # Deaths after stopping ART are HIV deaths
        return

    def make_modifiers(self, sim, hiv):
        """ ART's effects on HIV """
        st, a, p, name = self.graph.states, self.graph.attrs, self.pars, self.name
        on, ever = st.on_art, ~st.naive
        supp = a.suppressed
        eff = p.effective_art_efficacy*supp + p.nonsupp_art_efficacy*(1 - supp)
        v = self.value
        rel_eff, rel_m, rel_f = v('rel_art_mortality_effective'), v('rel_art_mortality_unsupp_m'), v('rel_art_mortality_unsupp_f')
        def art_mort(sim, art=self): # UPSTREAM: an amount mixing another module's attribute with ss.by() is evaluated for different agents (shapes differ), so it is a function here
            unsupp = np.where(sim.people.male.raw, rel_m, rel_f)
            return np.where(art.suppressed.raw > 0.5, rel_eff, unsupp)
        mods = [
            ss.modify(hiv.rel_trans, multiply=1 - eff*np.minimum(on.time/p.time_to_art_efficacy, 1), where=on, name=f'{name}_rel_trans'),
            ss.modify(hiv.progression, multiply=0, where=ever, name=f'{name}_halts_progression'), # Never resumes: ART replaces the natural history
            ss.modify(hiv.aids_death, multiply=0, where=ever, name=f'{name}_halts_aids_death'),
            set_attr(hiv.cd4, set=a.cd4_preart + (a.cd4_potential - a.cd4_preart)*(2/(1 + np.exp(-p.art_cd4_growth*on.time/ss.months(1))) - 1), where=on, name=f'{name}_cd4_on'), # As cd4_on
            set_attr(hiv.cd4, set=np.maximum(1, a.cd4_stop*(1 - st.off_art.progress)), where=st.off_art, name=f'{name}_cd4_off'),
            ss.modify(hiv.cd4_death, multiply=art_mort, where=on, name=f'{name}_mortality'),
            ss.modify(hiv.prenatal, multiply=1 - p.pmtct_efficacy, where=ss.mother(on), name=f'{name}_pmtct_prenatal'),
            ss.modify(hiv.postnatal, multiply=1 - p.pmtct_efficacy, where=ss.mother(on, born=True), name=f'{name}_pmtct_postnatal'),
        ]
        if self.art_death_age:
            bands = self.art_death_age
            mods.append(ss.modify(hiv.cd4_death, multiply=lambda sim: interval_mult(bands, sim.people.age.raw), where=on, name=f'{name}_mortality_age'))
        if self.art_death_dur:
            bands = self.art_death_dur
            def by_dur(sim, art=self):
                days = (sim.ti - art.ti_on_art.raw)*sim.t.dt.days
                return interval_mult(bands, np.nan_to_num(days, nan=-1))
            mods.append(ss.modify(hiv.cd4_death, multiply=by_dur, where=on, name=f'{name}_mortality_dur'))
        if self.value('rel_dur_on_art') != 1:
            mods.append(ss.modify(self.stop, multiply=1/self.value('rel_dur_on_art'), name=f'{name}_rel_dur'))
        if self.dur_on_art_trend is not None:
            trend = self.dur_on_art_trend
            mods.append(ss.modify(self.stop, multiply=lambda sim: 1/np.interp(sim.t.now('year'), trend['years'], trend['vals']), name=f'{name}_dur_trend'))
        return mods

    def ready(self, sim):
        """ Who can start ART now: diagnosed, infected, referred and past the delay, and naive or off ART """
        hiv = sim.diseases[self.disease]
        ti = sim.ti
        return hiv.diagnosed & hiv.infected & (self.ti_ready <= ti) & (self.naive | self.off_art)

    def ready_naive(self, sim):
        """ Who can start ART for the first time (with no coverage target, people who stop do not restart, as in STIsim 3) """
        return self.ready(sim) & self.naive

    def priority(self, sim, uids):
        """
        Who starts first under a coverage target: those who became ready on this step (just diagnosed, or past the
        referral delay) and take it up (art_initiation), then everyone else ready, lowest CD4 over care seeking first
        (care seeking is higher in pregnancy). As in STIsim 3, the newly ready are linked at once while there are places.
        """
        hiv = sim.diseases[self.disease]
        care = np.maximum(self.care_seeking.raw[uids], 0.01)
        preg = sim.demographics.get('pregnancy') if hasattr(sim, 'demographics') else None
        if preg is not None:
            care = care*np.where(preg.pregnant.raw[uids], self.maternal_care_scale, 1.0)
        score = -hiv.cd4_values(uids)/care
        p = self.pars.art_initiation
        newly = (self.ti_ready.raw[uids] == sim.ti) & (self.initiation.raw[uids] < float(getattr(p, 'raw', p)))
        return score + np.where(newly, 1e9, 0.0)

    def rank(self, sim, uids):
        return self.priority(sim, uids)

    def infected_pool(self, sim):
        """ People with HIV, whom a share of coverage is a share of """
        return sim.diseases[self.disease].infected

    def make_steps(self):
        """ Deliveries of ART, from the coverage: made here, so the sim has them from the start (see TargetProgram) """
        name = self.name
        tr = self.graph.transitions
        drug = ss.Treatment(moves=tr.start, name=f'{name}_drug')
        strata = coverage_keys(self.raw_coverage, ['n_art', 'p_art'])
        if not strata: # No target: everyone ready, once, with art_initiation
            p = self.pars.art_initiation
            return [ss.treat(drug, where=Lookup(name, 'ready_naive'), coverage=float(getattr(p, 'raw', p)), name=f'{name}_initiation')]
        out = []
        for i,stratum in enumerate(strata):
            among = Lookup(name, 'infected_pool', stratum)
            where = Lookup(name, 'ready', stratum)
            label = '' if stratum is None else f'_{i}'
            out.append(ss.treat(drug, where=where, target=stock_share(ShareLookup(name, i), among), priority=RankLookup(name),
                                above=tr.stop, name=f'{name}_delivery{label}'))
        return out

    def init_post(self):
        super().init_post()
        self.care_seeking.raw[:] = np.maximum(self.care_seeking.raw, 0.01) if hasattr(self.care_seeking, 'raw') else self.care_seeking
        hiv = self.sim.diseases[self.disease]
        dx = hiv.diagnosed.uids # Diagnosed at the start (init_diagnosed): ready at once
        self.ti_ready[dx[np.isnan(self.ti_ready[dx])]] = 0
        self.vls_coverage, self.vls_format, self.vls_age_bins, self.vls_sex_keys = parse_coverage(
            self.raw_vls_coverage, valid_names=['p_vls'], yearvec=self.t.yearvec, smoothness=self.smoothness, missing_fill=1.0)
        if self.vls_coverage is not None:
            vals = np.concatenate(list(self.vls_coverage.values())) if isinstance(self.vls_coverage, dict) else np.asarray(self.vls_coverage)
            if np.any(vals < 0) or np.any(vals > 1):
                errormsg = 'vls_coverage must be a proportion (0-1), not an absolute count.'
                raise ValueError(errormsg)
        return

    def step(self):
        self.correct_vls()
        ss.Disease.step(self)
        return

    def correct_vls(self):
        """ Keep the share suppressed among people on ART at vls_coverage, by adherence rank (stratified targets per stratum) """
        vls = self.vls_coverage
        if vls is None:
            return
        ti = self.ti
        on = self.on_art.uids
        if not len(on):
            return
        supp = self.suppressed.raw # UPSTREAM: a drawn attribute can only be changed by writing its array
        def fix(uids, p):
            if not len(uids):
                return
            n = int(np.clip(round(p*len(uids)), 0, len(uids)))
            order = uids[np.argsort(-self.adherence.raw[uids], kind='stable')]
            supp[order[:n]] = 1.0
            supp[order[n:]] = 0.0
        if isinstance(vls, dict):
            for ab in self.vls_age_bins:
                for sk in (self.vls_sex_keys or [None]):
                    key = (ab, sk) if sk is not None else ab
                    arr = vls.get(key, np.ones(len(self.t.yearvec)))
                    mask = np.asarray(age_sex_mask(ab, sk, self.sim.people).raw)[on]
                    fix(on[mask], arr[min(ti, len(arr) - 1)])
        else:
            fix(on, vls[min(ti, len(vls) - 1)])
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('p_on_art', dtype=float, label='Proportion of people with HIV on ART', scale=False),
            ss.Result('new_agents_on_art', dtype=float, label='New on ART', scale=False, auto_plot=False),
            ss.Result('n_on_art_pregnant', dtype=float, label='Pregnant women on ART', scale=False, auto_plot=False),
            ss.Result('n_vls', dtype=float, label='Number virally suppressed', scale=False, auto_plot=False),
            ss.Result('p_vls', dtype=float, label='Proportion of people with HIV virally suppressed', scale=False, auto_plot=False),
        )
        return

    def update_results(self):
        super().update_results()
        ti = self.ti
        res = self.results
        alive = self.sim.people.alive.raw
        on = self.on_art.raw & alive
        hiv = self.sim.diseases[self.disease]
        n_inf = np.count_nonzero(hiv.infected.raw & alive)
        n_vls = np.count_nonzero(on & (self.suppressed.raw > 0.5))
        w = ss.results.people_scale(self.sim) # A declared module's results are written in people
        res.p_on_art[ti] = sc.safedivide(np.count_nonzero(on), n_inf)
        res.n_vls[ti] = w*n_vls
        res.p_vls[ti] = sc.safedivide(n_vls, n_inf)
        preg = self.sim.demographics.get('pregnancy') if hasattr(self.sim, 'demographics') else None
        if preg is not None:
            res.n_on_art_pregnant[ti] = w*np.count_nonzero(on & preg.pregnant.raw)
        return

    def finalize_results(self):
        super().finalize_results()
        res = self.results
        res.new_agents_on_art.values[:] = res.start.values
        return


# %% Prevention products

class TargetProgram(ss.cascade):
    """
    A product delivered to a coverage target among the eligible, ranked by a willingness drawn for each person, with
    one delivery per stratum for stratified coverage; the base of `VMMC` and `Prep`

    The deliveries are made here, so the sim has them (and the products) from the start; what they read is
    looked up in the sim by name, since the sim works on copies.
    """
    valid_names = []
    stratum_sex = None # For age-only strata, the sex they are of

    def __init__(self, coverage=None, eligibility=None, smoothness=0, format_priority='n', name=None, label=None, disease='hiv'):
        self.raw_coverage = coverage
        self.eligibility = eligibility
        self.smoothness = smoothness
        self.format_priority = format_priority
        self.disease = disease
        self.prog_name = name if name is not None else type(self).__name__.lower()
        steps = self.make_steps()
        super().__init__(*steps, name=self.prog_name, label=label)
        self.eligibility = eligibility # ss.Intervention resets it
        self.define_states(ss.FloatArr('willingness', default=ss.random()))
        return

    def strata(self):
        """ The strata of the coverage: [], [None], or a list of (age bin, sex) """
        return coverage_keys(self.raw_coverage, self.valid_names, sex=self.stratum_sex)

    def make_steps(self):
        return []

    def target_deliveries(self, product, pool_name, where_name):
        """ One delivery per stratum, to the coverage target among the pool (the method named pool_name, of the sim) """
        out = []
        name = self.prog_name
        for i,stratum in enumerate(self.strata()):
            among = Lookup(name, pool_name, stratum)
            where = Lookup(name, where_name, stratum, also=pool_name)
            share = ShareLookup(name, i)
            label = '' if stratum is None else f'_{i}'
            out.append(ss.vaccinate(product, where=where, target=stock_share(share, among), priority=RankLookup(name),
                                    disease=self.disease, name=f'{name}_delivery{label}'))
        return out

    def init_pre(self, sim):
        super().init_pre(sim)
        strata, parsed = coverage_strata(self.raw_coverage, self.valid_names, sim, smoothness=self.smoothness,
                                         format_priority=self.format_priority, sex=self.stratum_sex)
        self.shares = [share for stratum,share in strata]
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parsed
        return

    def rank(self, sim, uids):
        return self.willingness.raw[uids]

    def held(self, product):
        """ Who holds a product, as a boolean array over all agents """
        ppl = self.sim.people
        arr = getattr(ppl, product.dim, None)
        if arr is None:
            return np.zeros(len(ppl.alive.raw), dtype=bool)
        return arr.raw > 0

    def step(self):
        return # The deliveries do the work


class Lookup:
    """ A method of a module in the sim, found by name when called (so it is the sim's copy), optionally restricted to an (age bin, sex) stratum and joined with another method """
    def __init__(self, modname, method, stratum=None, also=None):
        self.modname, self.method, self.stratum, self.also = modname, method, stratum, also
        return

    def __call__(self, sim):
        mod = find_mod(sim, self.modname)
        out = raw(getattr(mod, self.method)(sim))
        if self.also is not None:
            out = out & raw(getattr(mod, self.also)(sim))
        if self.stratum is not None:
            out = out & raw(age_sex_mask(*self.stratum, sim.people))
        return out


class RankLookup:
    """ The priority score of some agents from a module in the sim, found by name when called """
    def __init__(self, modname):
        self.modname = modname
        return

    def __call__(self, sim, uids):
        return find_mod(sim, self.modname).rank(sim, uids)


class ShareLookup:
    """ The share to cover in a stratum of a module's coverage, found by name when called """
    def __init__(self, modname, i):
        self.modname, self.i = modname, i
        return

    def __call__(self, sim, n_people):
        return find_mod(sim, self.modname).shares[self.i](sim, n_people)


class VMMC(TargetProgram):
    """
    Voluntary medical male circumcision, and traditional circumcision

    Program VMMC circumcises men up to a coverage target (a stock: circumcision prevalence among eligible men),
    highest willingness first; circumcision is never undone. Traditional circumcision is offered once to boys
    before sexual debut as they reach `traditional_age`, with probability `traditional_prob`. Circumcision
    reduces susceptibility to HIV by `eff_circ` (program) or `eff_circ_traditional`.

    Args:
        coverage: the coverage target (see `parse_coverage`; `n_vmmc`/`p_vmmc`, or stratified by age)
        eligibility (func): who program VMMC is for, a function of the sim (default: all men)
        eff_circ (float): reduction in susceptibility from program VMMC
        eff_circ_traditional (float): reduction in susceptibility from traditional circumcision
        traditional_prob (float): probability of traditional circumcision, as boys reach traditional_age
        traditional_age (float): the age of traditional circumcision

    Results: those of `ss.cascade`, and `new_circumcisions`, `n_circumcised`.

    Examples:
        ```python
        vmmc = sti.VMMC(coverage={'year': [2010, 2025], 'value': [0, 0.4]})
        both = sti.VMMC(coverage=0.3, eligibility=lambda sim: sim.people.age >= 25, traditional_prob=0.1)
        ```
    """
    valid_names = ['n_vmmc', 'p_vmmc']
    stratum_sex = 1

    def __init__(self, pars=None, coverage=None, eligibility=None, smoothness=0, format_priority='n', eff_circ=0.6,
                 eff_circ_traditional=0.6, traditional_prob=0, traditional_age=15, name=None, label=None, **kwargs):
        pars = sc.mergedicts(pars, kwargs)
        traditional_prob = pars.pop('traditional_prob', traditional_prob)
        self.traditional_prob = traditional_prob.pars.p if isinstance(traditional_prob, ss.bernoulli) else traditional_prob
        self.traditional_age = pars.pop('traditional_age', traditional_age)
        self.eff_circ = pars.pop('eff_circ', eff_circ)
        self.eff_circ_traditional = pars.pop('eff_circ_traditional', eff_circ_traditional)
        disease = pars.pop('disease', 'hiv')
        if pars:
            errormsg = f'Unknown parameters for VMMC: {sc.strjoin(pars.keys())}'
            raise ValueError(errormsg)
        super().__init__(coverage=coverage, eligibility=eligibility, smoothness=smoothness, format_priority=format_priority,
                         name=name, label=label, disease=disease)
        return

    def make_steps(self):
        name = self.prog_name
        self.program = ss.Vaccine(efficacy=self.eff_circ, name=f'{name}_program')
        self.traditional = ss.Vaccine(efficacy=self.eff_circ_traditional, name=f'{name}_traditional')
        steps = self.target_deliveries(self.program, 'pool', 'not_traditional')
        if self.traditional_prob:
            steps.append(ss.vaccinate(self.traditional, at_age=self.traditional_age, where=Lookup(name, 'before_debut'),
                                      coverage=self.traditional_prob, disease=self.disease, name=f'{name}_traditional_delivery'))
        return steps

    def pool(self, sim):
        """ Men eligible for program VMMC """
        out = sim.people.male & sim.people.alive
        return out & self.eligibility(sim) if self.eligibility is not None else out

    def not_traditional(self, sim):
        return ~self.held(self.traditional)

    def before_debut(self, sim):
        """ Boys before sexual debut, not circumcised by the program """
        net = sim.networks.get('structuredsexual')
        out = sim.people.male & ~self.held(self.program)
        return out & ~net.over_debut if net is not None else out

    @property
    def circumcised(self):
        """ Who is circumcised (by the program or traditionally), as a boolean array over all agents """
        return self.held(self.program) | self.held(self.traditional)

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_circumcisions', dtype=int, label='New circumcisions', auto_plot=False),
            ss.Result('n_circumcised', dtype=int, label='Number circumcised', auto_plot=False),
        )
        return

    def update_results(self):
        super().update_results()
        self.results.n_circumcised[self.ti] = np.count_nonzero(self.circumcised & self.sim.people.alive.raw)
        return

    def finalize_results(self):
        super().finalize_results()
        new = np.zeros(len(self.results.new_circumcisions))
        for step in self.steps:
            new += self.sim.results[step.name].delivered.values
        self.results.new_circumcisions.values[:] = new
        return


class Prep(TargetProgram):
    """
    Pre-exposure prophylaxis (PrEP)

    A PrEP product (oral or injectable) has an efficacy (`prep_eff`), an adherence (`prep_adh`, which multiplies
    the efficacy), and a course duration (`prep_dur`, the mean time until a course lapses). People are enrolled
    up to a coverage target (a stock: the share of the eligible and uninfected on PrEP), highest willingness
    first; someone on another PrEP product is not enrolled.

    Args:
        prep_eff (float): efficacy against acquisition, when fully adherent
        prep_dur (dur): mean duration of a course
        prep_adh (float): adherence (0-1), multiplied into the efficacy
        coverage: the coverage target (see `parse_coverage`; `n_prep`/`p_prep`, or stratified); default a ramp from 2004
        eligibility (func): who it is for, a function of the sim (default: female sex workers)

    Results: those of `ss.cascade`, and `new_prep`, `n_on_prep`.

    Examples:
        ```python
        prep = sti.Prep(coverage={'year': [2020, 2025], 'value': [0, 0.5]})
        prep = sti.Prep(coverage=0.4, eligibility=lambda sim: sim.people.female & (sim.people.age < 25))
        ```
    """
    valid_names = ['n_prep', 'p_prep']

    def __init__(self, pars=None, prep_eff=0.8, prep_dur=None, prep_adh=1.0, coverage=None, eligibility=None, smoothness=0,
                 format_priority='n', name=None, label=None, **kwargs):
        pars = sc.mergedicts(pars)
        years = kwargs.pop('years', None)
        if years is not None and coverage is not None:
            coverage = {'year': years, 'value': coverage}
        if coverage is None:
            coverage = {'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}
        if eligibility is None:
            eligibility = default_prep_eligibility
        self.prep_eff = pars.pop('prep_eff', prep_eff)
        self.prep_dur = pars.pop('prep_dur', prep_dur if prep_dur is not None else ss.months(3))
        self.prep_adh = pars.pop('prep_adh', prep_adh)
        super().__init__(coverage=coverage, eligibility=eligibility, smoothness=smoothness, format_priority=format_priority,
                         name=name, label=label, **kwargs)
        return

    def make_steps(self):
        self.product = ss.Vaccine(efficacy=self.prep_eff*self.prep_adh, waning=self.prep_dur, name=f'{self.prog_name}_product')
        return self.target_deliveries(self.product, 'pool', 'no_other')

    def pool(self, sim):
        """ The eligible and uninfected """
        hiv = sim.diseases[self.disease]
        return sim.people.alive & ~hiv.infected & self.eligibility(sim)

    def no_other(self, sim):
        """ Not on another PrEP product: one at a time """
        out = np.ones(len(sim.people.alive.raw), dtype=bool)
        for mod in sim.interventions.values():
            if isinstance(mod, Prep) and mod is not self:
                out &= ~mod.held(mod.product)
        return out

    @property
    def on_prep(self):
        """ Who is on this PrEP product, as a boolean array over all agents """
        return self.held(self.product)

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_prep', dtype=int, label='New PrEP starts', auto_plot=False),
            ss.Result('n_on_prep', dtype=int, label='Number on PrEP', auto_plot=False),
        )
        return

    def update_results(self):
        super().update_results()
        self.results.n_on_prep[self.ti] = np.count_nonzero(self.on_prep & self.sim.people.alive.raw)
        return

    def finalize_results(self):
        super().finalize_results()
        new = np.zeros(len(self.results.new_prep))
        for step in self.steps:
            new += self.sim.results[step.name].delivered.values
        self.results.new_prep.values[:] = new
        return


def default_prep_eligibility(sim):
    """ Female sex workers """
    return sim.networks.structuredsexual.fsw
