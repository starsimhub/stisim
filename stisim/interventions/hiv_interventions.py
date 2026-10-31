"""
Define HIV interventions for STIsim
By default, these all have units of a year and timesteps of 1/12
"""

import sciris as sc
import starsim as ss
import numpy as np
from stisim.diseases.hiv import art_pars
from stisim.interventions.base_interventions import STITest
from stisim.interventions.utils import (
    parse_coverage, compute_coverage_target, compute_stratum_targets,
    resolve_coverage_targets, age_sex_mask, top_n_by,
)
from stisim.utils import count


# %% HIV classes
__all__ = ["HIVDx", "HIVTest", "InfantHIVTest", "HIVCare", "ART", "VMMC", "Prep"]


class HIVDx(ss.Product):
    """HIV diagnostic product used by :class:`HIVTest`.

    A simple perfect-sensitivity test that classifies agents as positive
    (infected) or negative (susceptible) based on their current HIV state.

    Args:
        *args: Positional arguments forwarded to ``ss.Product``.
        **kwargs: Keyword arguments forwarded to ``ss.Product``.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.result_list = ['positive', 'negative']

    def administer(self, sim, uids):
        outcomes = {r: ss.uids() for r in self.result_list}
        outcomes['positive'] = sim.diseases.hiv.infected.uids.intersect(uids)
        outcomes['negative'] = sim.diseases.hiv.susceptible.uids.intersect(uids)
        return outcomes


class HIVTest(STITest):
    """
    HIV-specific testing intervention.

    Tests eligible agents for HIV; positive results set hiv.diagnosed=True,
    which is a prerequisite for ART initiation. By default, only undiagnosed
    agents are eligible.

    The testing → diagnosis → ART pipeline works as follows:

        1. HIVTest tests eligible agents each timestep (annual probability, converted via ss.probperyear)
        2. Positive results set hiv.diagnosed=True and hiv.ti_diagnosed
        3. HIVTest samples dur_dx2tx and refers them to treatment, to start at ti + delay (hiv.refer())
        4. ART picks up agents whose ti_art == current ti (and not yet on_art),
           filters by art_initiation, and puts them on ART
        5. If coverage data is provided, ART corrects to match targets

    Args:
        product: diagnostic product (default: :class:`HIVDx`).
        pars (dict): override default pars (``rel_test``, ``dt_scale``, ``dur_dx2tx``).
        test_prob_data: annual testing probability (if ``dt_scale=True``, the default).
            A value of 0.1 means ~10% of eligible agents tested per year. To
            specify a per-timestep probability instead, set ``dt_scale=False``.
        years (array): years over which testing is active (mutually exclusive with ``start``).
        start (float): calendar year when testing begins.
        eligibility (func): who can be tested. Default: undiagnosed agents.
        name, label: standard module identifiers.
        newborn_test (:class:`InfantHIVTest`): if supplied, schedules an infant
            HIV test at delivery for unborn children of mothers who test
            positive (requires ``ss.MaternalNet`` and ``ss.Pregnancy`` in the sim).
        dur_dx2tx (ss.Dist): delay from diagnosis to scheduled ART start
            (default: ``ss.constant(0)``, no delay). May be passed directly
            or via the ``pars`` dict. Example: ``ss.constant(ss.years(0.5))``
            for a 6-month delay.
        dt_scale (bool): if ``True`` (default), ``test_prob_data`` is an annual probability.
            Set to ``False`` for per-timestep probability.

    Example::

        # Test 20% of undiagnosed agents per year starting in 2000
        test = sti.HIVTest(test_prob_data=0.2, start=2000, name='hiv_test')

        # Test everyone every timestep (per-timestep probability)
        test = sti.HIVTest(test_prob_data=1.0, dt_scale=False, name='hiv_test')

        # FSW-targeted testing at higher rate
        fsw_test = sti.HIVTest(
            test_prob_data=0.5,
            name='fsw_test',
            eligibility=lambda sim: sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed,
        )

        # ANC testing: test undiagnosed pregnant women in first trimester
        anc_test = sti.HIVTest(
            test_prob_data=0.9,
            dt_scale=False,
            name='anc_test',
            eligibility=lambda sim: (sim.demographics.pregnancy.P & (sim.demographics.pregnancy.P.time < ss.weeks(13))
                                     & ~sim.diseases.hiv.diagnosed),
        )
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None,
                 eligibility=None, name=None, label=None, newborn_test=None, **kwargs):
        if product is None:
            suffix = f'_{name}' if name else ''
            product = HIVDx(name=f'HIVDx{suffix}')

        # Super init
        super().__init__(product=product, test_prob_data=test_prob_data, years=years,
                         start=start, eligibility=eligibility, name=name, label=label)

        # Deal with eligibility - default is undiagnosed people
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed

        # Newborn test
        self.newborn_test = newborn_test
        
        # Deal with pars
        self.define_pars(
            dur_dx2tx=ss.constant(0),
        )
        self.update_pars(pars, **kwargs)

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        hiv = sim.diseases.hiv
        hiv.diagnosed[pos_uids] = True
        hiv.ti_diagnosed[pos_uids] = self.ti
        # Schedule ART start for those not already on ART
        to_schedule = pos_uids[~hiv.on_art[pos_uids]]
        if len(to_schedule):
            delay = self.pars.dur_dx2tx.rvs(to_schedule, round=True)
            hiv.refer(to_schedule, self.ti + delay)

        # Schedule infant HIV test for unborn children of newly diagnosed mothers
        if self.newborn_test is not None and len(pos_uids):
            mn        = sim.get_module(ss.MaternalNet, die=False)
            pregnancy = sim.get_module(ss.Pregnancy, die=False)
            if mn is None or pregnancy is None:
                raise RuntimeError(
                    'HIVTest: cannot schedule newborn_test without ss.MaternalNet '
                    'and ss.Pregnancy in the sim'
                )
            pos_mask  = np.isin(mn.p1, pos_uids)
            unborn    = mn.p2[pos_mask]
            ti_births = np.ceil(pregnancy.ti_exit_P[mn.p1[pos_mask]]) # The step of the delivery
            valid     = ~np.isnan(ti_births)
            if valid.any():
                self.newborn_test.schedule(unborn[valid], ti_births[valid].astype(int))
        return outcomes


class InfantHIVTest(HIVTest):
    """
    HIV test for infants born to mothers diagnosed during pregnancy.

    Scheduled by :class:`HIVTest` (via ``newborn_test=``) or :class:`ANCTest`
    when the mother tests positive. Fires only at the scheduled timestep.
    A positive result sets ``hiv.diagnosed`` on the infant, enabling ART linkage.

    Args:
        test_prob (float): probability of the infant receiving the test once
                           scheduled (default 1.0 — if scheduled, it happens)
        name, label:       standard
    """
    def __init__(self, test_prob=1.0, name=None, label=None, **kwargs):
        super().__init__(name=name, label=label, **kwargs)
        self.test_prob_dist = ss.bernoulli(p=test_prob)

    def check_eligibility(self):
        # Only test infants at their scheduled timestep
        scheduled_now = (self.ti_scheduled == self.ti).uids
        return self.test_prob_dist.filter(scheduled_now)


class HIVCare(ss.Intervention):
    """
    HIV care: being in care after referral (e.g. a diagnosis), being lost to follow-up, and returning. With it, ART
    starts only for people in care (see sti.ART, which runs it at the start of its own step).

    Its defaults reproduce sti.ART without it: people referred (by sti.HIVTest, after its dur_dx2tx) are in care at once,
    and stay in care, through ART and after it.

    Each rate can be a rate (e.g. ``ss.peryear(0.1)``) or a function of the sim and the agents giving a rate per year for
    each (e.g. higher with AIDS).

    Args:
        dur_link (ss.Dist): time from referral to being in care, when ART can start
        rate_loss: loss to follow-up of people in care and not on ART
        rate_loss_art: loss to follow-up of people on ART, who stop it
        rate_return: return to care of people lost to follow-up, who can start ART again at once
        p_lost_at_stop (ss.bernoulli): of being lost to follow-up on stopping ART for any other reason (otherwise they stay in care, and can restart it)
    """
    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.define_pars(
            dur_link = ss.constant(0),
            rate_loss = ss.peryear(0),
            rate_loss_art = ss.peryear(0),
            rate_return = ss.peryear(0),
            p_lost_at_stop = ss.bernoulli(p=0),
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('in_care'),
            ss.BoolState('lost'), # Lost to follow-up
            ss.FloatArr('ti_link'), # When those referred will be in care
            ss.FloatArr('ti_in_care'),
            ss.FloatArr('ti_lost'),
        )
        self.u_loss = ss.random(name='u_loss')
        self.u_loss_art = ss.random(name='u_loss_art')
        self.u_return = ss.random(name='u_return')
        return

    def prob(self, key, uids):
        """ The probability in a step of HIV's of each agent's rate `key` """
        rate, dt = self.pars[key], self.sim.diseases.hiv.dt
        if callable(rate) and not isinstance(rate, ss.Rate):
            return -np.expm1(-np.asarray(rate(self.sim, uids), dtype=float)*dt.years)
        return rate.to_prob(dt)

    def refer(self, uids, ti, art):
        """ People referred at `ti` are in care after dur_link, when ART can start (ART's ti_art) """
        self.ti_link[uids] = ti + self.pars.dur_link.rvs(uids, round=True)
        art.ti_art[uids] = self.ti_link[uids]
        return

    def update(self, art):
        """ Linkage, loss to follow-up, and return to care, in this step (called by sti.ART) """
        ti = art.hiv.ti
        linked = (self.ti_link <= ti).uids
        self.in_care[linked] = True
        self.lost[linked] = False
        self.ti_in_care[linked] = ti
        self.ti_link[linked] = np.nan

        waiting = (self.in_care & ~art.on_treatment).uids
        lost = waiting[self.u_loss.rvs(waiting) < self.prob('rate_loss', waiting)]
        on_art = art.on_treatment.uids
        lost_art = on_art[self.u_loss_art.rvs(on_art) < self.prob('rate_loss_art', on_art)]
        if len(lost_art): # Only if anyone is: a draw of ART's in this step would change those after it
            art.stop_art(lost_art)
        self.lose(lost)
        self.lose(lost_art)

        lost_now = self.lost.uids
        back = lost_now[self.u_return.rvs(lost_now) < self.prob('rate_return', lost_now)]
        self.in_care[back] = True
        self.lost[back] = False
        self.ti_in_care[back] = ti
        art.ti_art[back] = ti # Ready to start ART again
        return

    def lose(self, uids):
        """ Lose people to follow-up """
        self.in_care[uids] = False
        self.lost[uids] = True
        self.ti_lost[uids] = self.sim.diseases.hiv.ti
        return

    def on_stop(self, uids):
        """ People stopping ART are lost to follow-up with probability p_lost_at_stop """
        self.lose(self.pars.p_lost_at_stop.filter(uids[self.in_care[uids]]))
        return

    def step(self):
        """ Run by sti.ART at the start of its step (see update()) """
        return

    def step_die(self, uids):
        self.in_care[uids] = False
        self.lost[uids] = False
        return


class ART(ss.Intervention):
    """
    Antiretroviral therapy intervention.

    Requires HIVTest (or equivalent) to diagnose agents first — ART only
    initiates agents who have hiv.diagnosed=True. A warning is raised if no
    HIVTest is found in the sim.

    Processing flow each timestep:
        1. Agents scheduled to stop ART are removed
        2. Agents whose scheduled ART start has arrived (ti_art == this timestep
           and not on_art) are filtered by art_initiation probability
        3. If coverage is specified: agents are added/removed to match the target
           number, prioritized by CD4 count and care-seeking propensity
        4. If no coverage is specified: all who pass art_initiation go on ART
           (no capacity constraint)
        5. Mothers on ART reduce infant susceptibility (prenatal via MaternalNet,
           postnatal via BreastfeedingNet) by pmtct_efficacy

    The scheduling (setting ti_art) is owned by the testing intervention that
    diagnoses the agent (HIVTest, ANCTest). ART here is passive: it picks up
    whoever is ready to start.

    Coverage is parsed by :func:`parse_coverage` and accepts:
        - ``None``: no coverage target; treat all who initiate (default)
        - Scalar (e.g. ``0.8``): constant proportion of infected on ART
        - Dict: ``{'year': [...], 'value': [...]}`` — interpolated to yearvec
        - Dict with mixed format: ``{'year': [...], 'value': [...], 'format': ['n','n','p','p']}``
        - Single-column DataFrame: index=years, column ``n_art`` or ``p_art``
        - Dual-column DataFrame: both ``n_art`` and ``p_art`` columns — uses
          ``format_priority`` to resolve (default: ``n_art`` when non-NaN)
        - Stratified DataFrame: columns Year, AgeBin (e.g. ``[15,25)``), and
          optionally Gender/Sex, plus a numeric value column

    Intervention ordering: HIVTest must appear before ART in the interventions
    list so that agents diagnosed this timestep can initiate ART in the same step.

    Args:
        coverage:         coverage target in any format above (default None)
        art_initiation:   probability a newly diagnosed person initiates ART
                          (default: ss.bernoulli(p=0.9)). Set to 1 to treat all
                          diagnosed.
        pmtct_efficacy:   efficacy of maternal ART in reducing infant susceptibility
                          to HIV (default 0.96). Applied to both prenatal (MaternalNet)
                          and postnatal (BreastfeedingNet) transmission. Set to 1.0 for
                          complete protection (previous default behavior).
        vls_coverage:     stock target for the fraction of agents on ART who are
                          virally suppressed. Accepts the same formats as
                          ``coverage`` (scalar, time-varying dict, DataFrame, or
                          age/sex-stratified DataFrame) but must be a proportion.
                          Corrected each step by :meth:`vls_stock_correction`.
                          Operates in tandem with ``p_effective_art`` (see
                          ``pars`` above), the same way ``coverage`` operates
                          alongside ``art_initiation``: ``p_effective_art`` sets
                          the per-initiate probability of starting suppressed,
                          ``vls_coverage`` reconciles the resulting stock to
                          observed data. Default ``None`` means no correction;
                          missing strata in a stratified input default to 100%.
        smoothness:       interpolation smoothness (0=linear, default)
        format_priority:  when both n_art and p_art are non-NaN, prefer this format
                          ('n' or 'p', default 'n')

    Examples::

        # Simple: 80% of infected on ART
        art = sti.ART(coverage=0.8)

        # 70% of newly-initiated agents achieve viral suppression
        art = sti.ART(coverage=0.8, vls_coverage=0.7)

        # Time-varying VLS coverage
        art = sti.ART(coverage=0.8, vls_coverage={'year': [2000, 2010, 2025], 'value': [0.5, 0.7, 0.9]})

        # Age/sex-stratified VLS coverage — e.g. lower suppression in young men
        vls_df = pd.DataFrame({
            'Year':    [2020, 2020, 2020, 2020],
            'AgeBin':  ['[15,25)', '[15,25)', '[25,100)', '[25,100)'],
            'Gender':  ['m', 'f', 'm', 'f'],
            'p_vls':   [0.5, 0.7, 0.75, 0.85],
        })
        art = sti.ART(coverage=0.8, vls_coverage=vls_df)

        # Time-varying coverage
        art = sti.ART(coverage={'year': [2000, 2010, 2025], 'value': [0, 0.5, 0.9]})

        # No coverage target — 90% of newly diagnosed initiate (default art_initiation)
        art = sti.ART()

        # From a CSV file
        art = sti.ART(coverage=pd.read_csv('art_coverage.csv').set_index('year'))

        # Historical n_art then projected p_art
        df = pd.read_csv('n_art.csv').set_index('year')
        df['p_art'] = np.nan
        df.loc[2023:, 'p_art'] = 0.90
        art = sti.ART(coverage=df)
    """

    def __init__(self, pars=None, coverage=None, vls_coverage=None, smoothness=0, format_priority='n', **kwargs):
        super().__init__()

        self.define_pars(
            art_initiation=ss.bernoulli(p=0.9),
            pmtct_efficacy=0.96,  # How much maternal ART reduces infant susceptibility
                                   # to HIV via MaternalNet (prenatal) and BreastfeedingNet
                                   # (postnatal). Conceptually like infant PrEP.
            p_effective_art=ss.bernoulli(p=1.0),  # Probability a newly-initiated agent achieves viral suppression

            # Infectiousness on ART
            effective_art_efficacy = 0.99, # Reduction in transmission on effective (virally suppressive) ART: PARTNER/PARTNER2/Opposites Attract ("U=U"), less a small residual
            nonsupp_art_efficacy = 0.35, # Reduction in transmission on non-suppressive ART: Quinn (Rakai 2000), applied to an assumed non-suppressed viral load distribution
            time_to_art_efficacy = ss.months(6), # Time to reach full ART efficacy: a linear increase in efficacy

            # Time on ART, and after it
            rate_art_dropout = ss.peryear(0), # Rate of stopping ART early (e.g. from poor adherence), besides at the end of dur_on_art; none by default
            dur_on_art = ss.lognorm_ex(ss.years(3), ss.years(1.5)), # Base ART duration (scaled by rel_dur_on_art and trend)
            rel_dur_on_art = 1.0, # Calibratable scalar that scales ART duration
            dur_on_art_trend = None, # Optional time-varying scale, e.g. sc.objdict(years=np.array([2004, 2015, 2030]), vals=np.array([0.5, 1.0, 1.5]))
            dur_post_art = ss.normal(), # Time from stopping ART until the CD4 count reaches zero; set for each person as they stop (see stop_art()), in years
            dur_post_art_scale_factor = 0.1,

            # CD4 count on ART
            art_cd4_growth = 0.1, # Unitless parameter of how CD4 reconstitutes after starting ART, in a logistic growth function
            art_cd4_pars = dict(cd4_max=1000, cd4_healthy=500),

            # Mortality on ART, relative to the off-ART CD4-based rate, so that being on ART is never more dangerous than being
            # off it at the same CD4 count (as long as the largest of these times the largest age and duration multipliers is
            # at most 1, which init_pre() checks)
            rel_art_mortality_effective = 0.25, # Fraction of off-ART mortality retained on effective (suppressive) ART
            rel_art_mortality_unsupp_m = 0.7, # ... on non-suppressive ART, males
            rel_art_mortality_unsupp_f = 0.35, # ... females: the non-suppressive/effective ratio is 2x higher for men than for women
            art_death_age = [(0, 25, 1.0), (25, 35, 1.10), (35, 45, 1.21), (45, 125, 1.32)], # (age_lo, age_hi, mult)
            art_death_dur = None, # (dur_lo, dur_hi, mult) list, days since ti_art; None = no-op

            # Care seeking, which sets priority for ART
            care_seeking = ss.normal(loc=1, scale=0.5), # Relative care-seeking behavior (wide spread for low-care-seeking extremes)
            maternal_care_scale = 2, # Factor for scaling up care-seeking behavior during pregnancy
            art_priority = None, # Under a coverage target, a function of the sim and agents giving each a priority (higher starts first, and lower stops first); by default, lowest CD4 count over care seeking

            # Viral suppression on ART, after its start (each a rate, or a function of the sim and agents giving a rate per year for each)
            rate_failure = ss.peryear(0), # Of becoming unsuppressed
            rate_regain = ss.peryear(0), # Of becoming suppressed again
        )
        self.update_pars(pars, **kwargs)
        self.given_pars = set(sc.mergedicts(pars, kwargs)) # To tell those given here from those given to HIV (see init_pre())

        self._raw_coverage     = coverage
        self._raw_vls_coverage = vls_coverage
        self._smoothness       = smoothness
        self._format_priority  = format_priority
        self.coverage          = None  # Set in init_pre
        self.coverage_format   = None  # 'n', 'p', or per-timestep array
        self.age_bins          = None  # For stratified coverage
        self.sex_keys          = None  # For stratified coverage
        self.vls_coverage      = None  # Set in init_pre
        self.vls_format        = None  # 'p' or per-timestep array
        self.vls_age_bins      = None  # For stratified VLS coverage
        self.vls_sex_keys      = None  # For stratified VLS coverage

        # States
        self.define_states(
            ss.FloatArr('adherence', default=ss.random()),
            ss.BoolState('on_effective_art'), # On ART and virally suppressed
            ss.BoolState('on_nonsuppressive_art'), # On ART but not suppressed
            ss.BoolState('art_discontinued'), # Off ART, having been on it
            ss.BoolState('art_naive', default=True), # Never initiated ART
            ss.FloatArr('ti_art'), # Time of starting ART, or of a scheduled start (see sti.HIVTest)
            ss.TimeArr('ti_stop_art'),
            ss.FloatArr('cd4_preart'), # CD4 immediately before initiating ART
            ss.FloatArr('cd4_potential'), # Potential CD4 count if continually treated
            ss.FloatArr('cd4_postart'), # CD4 after stopping ART
            ss.FloatArr('dur_post_art'), # The time from stopping ART until the CD4 count reaches zero, in steps
            ss.FloatArr('baseline_care_seeking'),
            ss.FloatArr('care_seeking'),
        )
        self.u_dropout = ss.random(name='u_dropout')
        self.u_failure = ss.random(name='u_failure')
        self.u_regain = ss.random(name='u_regain')

        return

    def init_pre(self, sim):
        # ART's parameters given to HIV (as in STIsim 3) are taken from there, unless given here too
        hiv = sim.diseases.hiv
        from_hiv = {key: val for key, val in hiv.pars.items() if key in art_pars and val is not None}
        both = [key for key in from_hiv if key in self.given_pars]
        if both:
            errormsg = f'ART parameters given to both HIV and ART: {sc.strjoin(both)}. Give them to ART only.'
            raise ValueError(errormsg)
        for key, val in from_hiv.items():
            self.pars[key] = sc.dcp(val)
            hiv.pars[key] = None # Now ART's
        super().init_pre(sim)
        hiv.add_treatment(self)
        self.hiv = hiv
        cares = [m.name for m in sim.interventions() if isinstance(m, HIVCare)]
        if len(cares) > 1:
            raise ValueError(f'At most one sti.HIVCare, not {len(cares)}')
        self.care_name = cares[0] if cares else None

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

        # Warn if no HIV testing intervention found, or if it comes after ART
        intvs = list(sim.interventions())
        hiv_test_indices = [i for i, m in enumerate(intvs) if isinstance(m, HIVTest)]
        art_index = next((i for i, m in enumerate(intvs) if m is self), len(intvs))
        if not hiv_test_indices:
            ss.warn('ART intervention added without an HIV testing intervention; diagnosed agents may never be identified')
        elif all(i > art_index for i in hiv_test_indices):
            ss.warn('HIVTest appears after ART in the intervention list; agents diagnosed this timestep will not initiate ART until the next step')

        # Parse coverage data
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parse_coverage(
            self._raw_coverage, valid_names=['n_art', 'p_art'], yearvec=self.t.yearvec,
            smoothness=self._smoothness, format_priority=self._format_priority,
        )

        # missing_fill=1.0 so strata absent from a stratified vls_coverage default
        # to 100% suppressed (opposite of the 0% default used for `coverage` above).
        self.vls_coverage, self.vls_format, self.vls_age_bins, self.vls_sex_keys = parse_coverage(
            self._raw_vls_coverage, valid_names=['p_vls'], yearvec=self.t.yearvec,
            smoothness=self._smoothness, missing_fill=1.0,
        )
        if self.vls_coverage is not None:
            # Check actual values, not the format tag: parse_coverage always tags scalar
            # and single-column-DataFrame inputs 'p' regardless of magnitude, so an
            # out-of-range value (e.g. a percentage instead of a fraction) would slip
            # past a format=='n' check.
            vls_values = (np.concatenate(list(self.vls_coverage.values()))
                          if isinstance(self.vls_coverage, dict) else self.vls_coverage)
            if np.any(vls_values < 0) or np.any(vls_values > 1):
                errormsg = 'vls_coverage must be a proportion (0-1), not an absolute count.'
                raise ValueError(errormsg)

        self.initialized = True
        return

    # Treatment of HIV: the interface it calls (see sti.HIV.add_treatment())
    @property
    def on_treatment(self):
        """ On ART, suppressed or not """
        return self.on_effective_art | self.on_nonsuppressive_art

    @property
    def care(self):
        """ The care module (sti.HIVCare), if any: then ART starts only for those in care """
        return self.sim.interventions[self.care_name] if self.care_name else None

    @property
    def in_care(self):
        """ Who can start ART: those in care, if there is a care module, else everyone """
        care = self.care
        return care.in_care if care is not None else self.hiv.infected | ~self.hiv.infected

    def refer(self, uids, ti):
        """ Schedule ART to start at `ti` for those referred (e.g. diagnosed) and not on it already: through care, if there is a care module """
        uids = uids[~self.on_treatment[uids]]
        if self.care is not None:
            self.care.refer(uids, ti, self)
        else:
            self.ti_art[uids] = ti
        return

    def update_cd4(self):
        """ CD4 counts: rising on ART, and falling to 1 after it """
        hiv, ti = self.hiv, self.hiv.ti
        on_art = self.on_treatment.uids
        hiv.cd4[on_art] = self.cd4_increase(on_art)
        stopped = self.art_discontinued.uids
        start, end = self.ti_stop_art[stopped], hiv.ti_zero[stopped]
        share = np.clip(np.divide(ti - start, end - start, out=np.ones(len(start)), where=end > start), 0, 1) # Of the time to a CD4 count of zero elapsed
        hiv.cd4[stopped] = np.maximum(1, self.cd4_postart[stopped] - (self.cd4_postart[stopped] - 1)*share)
        return

    def step_treatment(self):
        """ Stopping ART early, at rate_art_dropout """
        on_art = self.on_treatment.uids
        rate = self.pars.rate_art_dropout
        self.stop_art(on_art[self.u_dropout.rvs(on_art) < rate.to_prob(self.hiv.dt)])
        return

    def death_multiplier(self, uids):
        """ The multiplier of ART on the death rates of some infected agents: rel_art_mortality (by suppression, and by sex if not suppressed), by age (art_death_age), and by time on ART (art_death_dur); 1 for those not on it """
        mult = np.ones(len(uids))
        on = self.on_treatment[uids]
        mult[on] = self.on_art_death_rate(uids[on])/np.maximum(self.hiv.off_art_death_rate(uids[on]), 1e-300)
        return mult

    def trans_multiplier(self, uids):
        """ The multiplier of ART on the infectiousness of some agents on it """
        return self.art_factor(uids)

    def step_die(self, uids):
        """ The people who died leave the states of ART, and of care """
        for state in [self.on_effective_art, self.on_nonsuppressive_art, self.art_discontinued]:
            state[uids] = False
        if self.care is not None:
            self.care.step_die(uids)
        return

    def init_care_seeking(self):
        """ Set the care-seeking behavior of new agents """
        uids = ss.uids(self.care_seeking.isnan)
        self.care_seeking[uids] = self.pars.care_seeking.rvs(uids)
        self.care_seeking[uids] = np.maximum(self.care_seeking[uids], 0.01) # Floor to prevent division by zero
        self.baseline_care_seeking[uids] = sc.dcp(self.care_seeking[uids]) # Copy it so pregnancy can modify it
        return

    def update_care_seeking(self):
        """ Care seeking of new agents; pregnant women seek care more, which makes them much less likely to stop treatment """
        self.init_care_seeking()
        if self.hiv.include_mtct:
            pregnant = self.sim.demographics.pregnancy.P
            self.care_seeking[pregnant] = self.baseline_care_seeking[pregnant] * self.pars.maternal_care_scale
            self.care_seeking[~pregnant] = self.baseline_care_seeking[~pregnant]
        return

    def init_post(self):
        super().init_post()
        self.init_care_seeking()
        return

    def on_art_death_rate(self, uids):
        """
        The annual rate of death of each agent on ART: the off-ART rate at the same CD4 count, times rel_art_mortality
        (by suppression, and by sex if not suppressed) and the multipliers by age (art_death_age) and time on ART (art_death_dur)
        """
        hiv = self.hiv
        male = self.sim.people.male[uids]
        effective = self.on_effective_art[uids]
        age = self.sim.people.age[uids]
        rate = hiv.off_art_death_rate(uids)
        rate[effective] *= self.pars.rel_art_mortality_effective
        rate[~effective & male] *= self.pars.rel_art_mortality_unsupp_m
        rate[~effective & ~male] *= self.pars.rel_art_mortality_unsupp_f

        age_mult = np.ones(len(uids))
        for age_lo, age_hi, mult in self.pars.art_death_age:
            age_mult[(age >= age_lo) & (age < age_hi)] *= mult
        rate = rate * age_mult

        if self.pars.art_death_dur is not None:
            dur_days = (hiv.ti - self.ti_art[uids]) * hiv.dt.days
            dur_mult = np.ones(len(uids))
            for dur_lo, dur_hi, mult in self.pars.art_death_dur:
                dur_mult[(dur_days >= dur_lo) & (dur_days < dur_hi)] *= mult
            rate = rate * dur_mult
        return rate

    def art_factor(self, uids):
        """ The multiplier of ART on the infectiousness of some agents on it, which reaches its full efficacy by a linear rise over time_to_art_efficacy """
        p, hiv = self.pars, self.hiv
        full_eff = np.where(self.on_effective_art[uids], p.effective_art_efficacy, p.nonsupp_art_efficacy)
        return 1 - full_eff*np.minimum((hiv.ti - self.ti_art[uids])/(p.time_to_art_efficacy/hiv.dt), 1)

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
        dur_art = self.hiv.ti - self.ti_art[uids]
        cd4_preart = self.cd4_preart[uids]
        cd4_total_gain = self.cd4_potential[uids] - cd4_preart
        return 2*cd4_total_gain/(1 + np.exp(-dur_art*self.pars.art_cd4_growth)) - cd4_total_gain + cd4_preart # Concave logistic

    def start_art(self, uids, p_effective_art=None):
        """
        Put agents on ART: suppressed, or not; it takes them out of HIV's natural history

        Args:
            uids: agents starting ART treatment
            p_effective_art (float/ss.Dist, optional): probability that each agent achieves viral suppression (effective
                ART) rather than non-suppressive ART. Pass a float to override `self.pars.p_effective_art`'s probability,
                or an already-initialized `ss.Dist` to use it directly. Defaults to `self.pars.p_effective_art` (default:
                always effective).
        """
        hiv = self.hiv
        ti = hiv.ti
        p = self.pars
        newly_treated = uids[self.art_naive[uids]]
        self.art_naive[newly_treated] = False
        self.ti_art[uids] = ti
        self.cd4_preart[uids] = hiv.cd4[uids]

        # Whether each achieves viral suppression (effective ART) or not, e.g. due to poor adherence or drug resistance
        if isinstance(p_effective_art, ss.Dist):
            effective = p_effective_art.filter(uids)
        else:
            if p_effective_art is not None:
                p.p_effective_art.set(p=p_effective_art)
            effective = p.p_effective_art.filter(uids)
        for state in [hiv.acute, hiv.latent, hiv.falling, self.art_discontinued]:
            state[uids] = False
        self.on_effective_art[uids] = False
        self.on_effective_art[effective] = True
        self.on_nonsuppressive_art[uids] = True
        self.on_nonsuppressive_art[effective] = False
        later = uids[hiv.ti_zero[uids] > ti] # Those due to die of AIDS now still do
        hiv.ti_zero[later] = np.nan
        hiv.ti_latent[uids] = hiv.ti_falling[uids] = np.nan

        # When they will stop
        dur_on_art = p.dur_on_art.rvs(uids) * p.rel_dur_on_art
        trend = p.dur_on_art_trend
        if trend is not None:
            dur_on_art = dur_on_art * np.interp(hiv.t.now('year'), trend.years, trend.vals)
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
        hiv = self.hiv
        ti = hiv.ti
        uids = (self.on_treatment & self.ti_stop_art.due()).uids if uids is None else ss.uids(uids)
        on = uids[self.on_treatment[uids]]
        factor = self.art_factor(on)
        hiv.rel_trans[on] = np.divide(hiv.rel_trans[on], factor, out=np.ones(len(on)), where=factor > 0) # Off ART from this step's transmission on, as they are put on it (see sti.HIV.update_transmission())
        self.on_effective_art[uids] = False
        self.on_nonsuppressive_art[uids] = False
        self.art_discontinued[uids] = True
        self.ti_stop_art[uids] = ti
        self.cd4_postart[uids] = hiv.cd4[uids]
        dur_mean = np.log(self.cd4_preart[uids])*hiv.cd4[uids]/self.cd4_potential[uids]
        dur_scale = dur_mean*self.pars.dur_post_art_scale_factor
        self.pars.dur_post_art.set(loc=ss.years(dur_mean), scale=ss.years(np.maximum(dur_scale, 1e-3))) # Not zero
        dur_post_art = self.pars.dur_post_art.rvs(uids)
        if np.isnan(dur_post_art).any():
            bad = uids[np.isnan(dur_post_art)]
            errormsg = f'Invalid post-ART duration for {bad}: CD4 {hiv.cd4[bad]}, before ART {self.cd4_preart[bad]}, potential {self.cd4_potential[bad]}'
            raise ValueError(errormsg)
        self.dur_post_art[uids] = np.maximum(dur_post_art.astype(int), 0)
        hiv.ti_zero[uids] = ti + self.dur_post_art[uids]
        if self.care is not None:
            self.care.on_stop(uids)
        return

    def set_suppression(self, uids, effective):
        """ Make agents on ART virally suppressed (effective=True), or not """
        uids = ss.uids(uids)
        uids = uids[self.on_treatment[uids]]
        self.on_effective_art[uids] = effective
        self.on_nonsuppressive_art[uids] = not effective
        return

    def _suppress_to_target(self, hiv, pool, target):
        """
        Suppress the top ``target`` agents in ``pool`` by ``adherence``;
        un-suppress the rest.
        """
        uids = pool.uids
        if not len(uids):
            return
        n_target = int(np.clip(round(target), 0, len(uids)))
        order = uids[np.argsort(-self.adherence[uids])]
        keep, drop = order[:n_target], order[n_target:]
        self.set_suppression(keep, effective=True)
        self.set_suppression(drop, effective=False)

    def vls_stock_correction(self, sim):
        """
        Correct viral suppression among agents already on ART toward the
        ``vls_coverage`` target. Stratified targets are corrected per stratum
        so age/sex differentials in the input survive the ranking.
        """
        if self.vls_coverage is None:
            return

        hiv = sim.diseases.hiv
        pool = hiv.on_art
        if not pool.any():
            return

        total, stratum_targets = resolve_coverage_targets(
            self.vls_coverage, self.vls_format, self.vls_age_bins, self.vls_sex_keys,
            self.ti, pool.uids, sim,
        )
        if stratum_targets is not None:
            for key, target in stratum_targets.items():
                ab, sex = key if isinstance(key, tuple) else (key, None)
                stratum = pool & age_sex_mask(ab, sex, sim.people)
                self._suppress_to_target(hiv, stratum, target)
        elif total is not None:
            self._suppress_to_target(hiv, pool, total)

    def step(self):
        """
        Apply ART at each timestep: stop ART for those scheduled, initiate for newly diagnosed,
        and correct coverage to match targets (per-stratum when stratified data is supplied,
        otherwise against the aggregate total).
        """
        sim = self.sim
        hiv = sim.diseases.hiv
        inf_uids = hiv.infected.uids
        self.update_care_seeking()
        if self.care is not None:
            self.care.update(self)
        self.step_suppression()

        # Determine treatment target (None = no capacity constraint). When
        # stratum_targets is set, correct within each stratum so age/sex
        # differentials aren't washed out by a global CD4-priority allocation.
        n_to_treat, stratum_targets = resolve_coverage_targets(
            self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
            self.ti, inf_uids, self.sim,
        )

        # Check who is stopping ART
        if hiv.on_art.any():
            stopping = hiv.on_art & (self.ti_stop_art <= self.ti)
            if stopping.any():
                try:
                    self.stop_art(stopping.uids)
                except:
                    errormsg = f'Error stopping ART for {stopping.uids}'
                    raise ValueError(errormsg)

        # Initiate ART for agents whose scheduled start time has arrived
        ready = ((self.ti_art == self.ti) & ~hiv.on_art & self.in_care).uids
        if len(ready):
            dx_to_treat = self.pars.art_initiation.filter(ready)

            if n_to_treat is None:
                # No coverage target — treat all who initiate
                self.start_art(dx_to_treat, p_effective_art=self.pars.p_effective_art)
            else:
                # Coverage target — only treat if spots available
                on_art = hiv.on_art
                n_available_spots = n_to_treat - len(on_art.uids)
                if n_available_spots > 0:
                    self.prioritize_art(sim, n=n_available_spots, awaiting_art_uids=dx_to_treat)

        # Correct coverage to match target (only if target is set).
        # Stratified targets get per-stratum correction; aggregate gets global.
        if n_to_treat is not None:
            self.art_coverage_correction(sim, target_coverage=n_to_treat,
                                         stratum_targets=stratum_targets)

        # Run after coverage correction so agents added to ART this step are in the pool.
        self.vls_stock_correction(sim)

        # PMTCT: reduce susceptibility of infants whose mothers are on ART.
        # This applies to both prenatal (MaternalNet) and postnatal (BreastfeedingNet)
        # transmission. Conceptually, maternal ART acts like infant PrEP — viral
        # suppression in the mother reduces the infant's exposure to HIV.
        # Note: the mother's own rel_trans is also reduced by ART (in update_transmission),
        # so total protection compounds both effects.
        pmtct_eff = self.pars.pmtct_efficacy
        if hasattr(sim.people, 'pregnancy'):
            preg = sim.people.pregnancy

            # Prenatal: protect unborn infants of pregnant mothers on ART
            if hasattr(sim.networks, 'maternalnet'):
                art_pregnant = (hiv.on_art & preg.P).uids
                if len(art_pregnant):
                    unborn = sim.networks.maternalnet.find_contacts(art_pregnant)
                    hiv.rel_sus[ss.uids(unborn)] *= 1 - pmtct_eff

            # Postnatal: protect breastfed infants of breastfeeding mothers on ART
            if hasattr(sim.networks, 'breastfeedingnet'):
                art_breastfeeding = (hiv.on_art & preg.lactating).uids
                if len(art_breastfeeding):
                    breastfed = sim.networks.breastfeedingnet.find_contacts(art_breastfeeding)
                    hiv.rel_sus[ss.uids(breastfed)] *= 1 - pmtct_eff

        return

    def priority(self, uids):
        """ Each agent's priority for ART under a coverage target (higher starts first, and lower stops first): art_priority, or by default the lowest CD4 count over care seeking """
        if self.pars.art_priority is not None:
            return np.asarray(self.pars.art_priority(self.sim, uids), dtype=float)
        return -self.hiv.cd4[uids]/self.care_seeking[uids]

    def step_suppression(self):
        """ Becoming unsuppressed (rate_failure), and suppressed again (rate_regain), on ART """
        for rate, dist, frm, effective in [('rate_failure', self.u_failure, self.on_effective_art, False), ('rate_regain', self.u_regain, self.on_nonsuppressive_art, True)]:
            uids = frm.uids
            par = self.pars[rate]
            prob = -np.expm1(-np.asarray(par(self.sim, uids), dtype=float)*self.hiv.dt.years) if callable(par) and not isinstance(par, ss.Rate) else par.to_prob(self.hiv.dt)
            if np.any(prob):
                self.set_suppression(uids[dist.rvs(uids) < prob], effective)
        return

    def prioritize_art(self, sim, n=None, awaiting_art_uids=None):
        """ Prioritize ART to n agents among those awaiting treatment """
        hiv = sim.diseases.hiv
        if awaiting_art_uids is None:
            # Agents whose scheduled ART start has arrived but who aren't on ART yet
            awaiting_art_uids = ((self.ti_art <= self.ti) & ~hiv.on_art & self.in_care).uids

        # Enough spots for everyone
        if n > len(awaiting_art_uids):
            start_uids = awaiting_art_uids

        # Not enough spots — prioritize lowest-CD4, highest-care-seeking
        else:
            start_uids = top_n_by(awaiting_art_uids, self.priority(awaiting_art_uids), n)

        self.start_art(start_uids, p_effective_art=self.pars.p_effective_art)

        return

    def art_coverage_correction(self, sim, target_coverage=None, stratum_targets=None):
        """
        Adjust ART coverage to match data.

        If ``stratum_targets`` is supplied, correction is done independently
        within each ``(age_bin, sex)`` stratum so that age/sex differentials
        in the input data are preserved. Otherwise, correction is done
        globally against ``target_coverage``.
        """
        if stratum_targets is not None:
            self._art_coverage_correction_stratified(sim, stratum_targets)
        else:
            self._art_coverage_correction_global(sim, target_coverage)
        return

    def _art_coverage_correction_global(self, sim, target_coverage):
        """ Adjust ART coverage globally to match an aggregate target """
        hiv = sim.diseases.hiv
        on_art = hiv.on_art

        # Too many on treatment → remove highest-CD4, lowest-care-seeking
        if len(on_art.uids) > target_coverage:
            n_to_stop    = int(len(on_art.uids) - target_coverage)
            on_art_uids  = on_art.uids
            stop_uids    = top_n_by(on_art_uids, -self.priority(on_art_uids), n_to_stop)
            self.ti_stop_art[stop_uids] = self.ti
            self.stop_art(stop_uids)

        # Not enough on treatment → add (only those whose scheduled start has arrived)
        elif len(on_art.uids) < target_coverage:
            n_to_add = target_coverage - len(on_art.uids)
            awaiting_art_uids = ((self.ti_art <= self.ti) & ~hiv.on_art & self.in_care).uids
            self.prioritize_art(sim, n=n_to_add, awaiting_art_uids=awaiting_art_uids)

        return

    def _art_coverage_correction_stratified(self, sim, stratum_targets):
        """
        Adjust ART coverage independently within each (age_bin, sex) stratum.

        This preserves the age/sex differentials in the stratified coverage
        input — e.g. young men and older women each reach their own target
        rather than being averaged out by a global CD4-priority allocation.
        Within a stratum, removals still drop the highest-CD4 agents first
        and additions still prioritize the lowest-CD4 agents (same priority
        rule as the global path).
        """
        hiv = sim.diseases.hiv
        ppl = sim.people

        for key, target_n in stratum_targets.items():
            # key is (age_bin, sex) when sex-stratified, else age_bin only
            if isinstance(key, tuple):
                ab, sex = key
            else:
                ab, sex = key, None

            in_stratum = age_sex_mask(ab, sex, ppl) & ppl.alive & hiv.infected
            on_art_in_stratum = (in_stratum & hiv.on_art).uids
            n_on_art = len(on_art_in_stratum)

            # Too many on treatment in this stratum → remove highest-CD4
            if n_on_art > target_n:
                n_to_stop    = n_on_art - target_n
                stop_uids    = top_n_by(on_art_in_stratum, -self.priority(on_art_in_stratum), n_to_stop)
                self.ti_stop_art[stop_uids] = self.ti
                self.stop_art(stop_uids)

            # Not enough in this stratum → add from diagnosed in same stratum
            elif n_on_art < target_n:
                n_to_add = target_n - n_on_art
                awaiting_art_uids = (in_stratum & hiv.diagnosed & ~hiv.on_art).uids
                self.prioritize_art(sim, n=n_to_add, awaiting_art_uids=awaiting_art_uids)

        return


class VMMC(ss.Intervention):
    """
    Voluntary medical male circumcision.

    Reduces male susceptibility to HIV acquisition (see HIV's ``eff_circ``).
    Unlike ART, VMMC does not require diagnosis — it circumcises males up to a
    coverage target, prioritized by willingness (a random per-agent score).

    If no coverage is specified, VMMC does nothing. Coverage must be provided
    explicitly via the coverage parameter.

    Coverage is parsed by :func:`parse_coverage` (same formats as ART, using
    ``n_vmmc``/``p_vmmc`` column names). Age-only stratification is supported
    (no Gender column required, since VMMC is males-only).

    Circumcision status (``circumcised``/``circ_traditional``/``ti_circumcised``) and
    efficacy (``eff_circ``) live on :class:`~stisim.diseases.hiv.HIV`, not on
    this intervention — see ``HIV.circumcise()`` and HIV's ``eff_circ``. This
    lets other circumcision sources (e.g. non-program/traditional circumcision)
    share the same state and be included in this intervention's "already
    circumcised" checks without going through VMMC at all.

    In addition to program VMMC (the coverage-target logic above), this class
    can optionally also model non-program/traditional circumcision.

    Args:
        coverage:          coverage target (default None; program VMMC does nothing
                           without data). See :func:`parse_coverage` for supported formats.
        eligibility:       optional function to restrict who is eligible for PROGRAM VMMC
                           (default: all males); does not affect traditional circumcision.
        smoothness:        interpolation smoothness (0=linear, default)
        format_priority:   when both n_vmmc and p_vmmc are non-NaN, prefer this format
        traditional_prob:  probability (float or ``ss.bernoulli``) that a pre-sexual-debut
                           male receives traditional circumcision once he reaches
                           ``traditional_age`` (default 0 = off).
        traditional_age:   age (years) at which pre-debut males are assessed for
                           traditional circumcision, once (default 15).

    Examples::

        # Program VMMC
        vmmc = sti.VMMC(coverage=0.3)
        vmmc = sti.VMMC(coverage={'year': [2010, 2025], 'value': [0, 0.4]})

        # Traditional circumcision only, no program VMMC
        trad = sti.VMMC(traditional_prob=0.1, traditional_age=15)

        # Program VMMC for adults + traditional circumcision for boys, combined
        both = sti.VMMC(coverage=0.3, eligibility=lambda sim: sim.people.age >= 25,
                        traditional_prob=0.1, traditional_age=15)
    """

    def __init__(self, pars=None, coverage=None, eligibility=None, smoothness=0, format_priority='n', **kwargs):
        super().__init__(eligibility=eligibility)

        self.define_pars(
            traditional_prob=ss.bernoulli(p=0),
            traditional_age=15,
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._smoothness      = smoothness
        self._format_priority = format_priority
        self.coverage         = None
        self.coverage_format  = None
        self.age_bins         = None
        self.sex_keys         = None

        # States
        self.define_states(
            ss.FloatArr('willingness', default=ss.random()),
            ss.BoolArr('traditional_assessed', default=False),
        )
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parse_coverage(
            self._raw_coverage, valid_names=['n_vmmc', 'p_vmmc'], yearvec=self.t.yearvec,
            smoothness=self._smoothness, format_priority=self._format_priority,
        )
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_circumcisions', dtype=int, label='New circumcisions', auto_plot=False),
            ss.Result('n_circumcised',     dtype=int, label='Number circumcised', auto_plot=False),
        )
        return

    def _circumcise_to_target(self, hiv, pool, target):
        """Top ``pool`` (a male BoolArr) up to ``target`` circumcised, choosing
        the highest-willingness uncircumcised men. Never removes (circumcision
        is irreversible). Returns the number of new circumcisions."""
        n_add = int(target) - (pool & hiv.circumcised).count()
        if n_add <= 0:
            return 0
        candidates = (pool & ~hiv.circumcised).uids
        if len(candidates) == 0:
            return 0
        new_circs = top_n_by(candidates, self.willingness[candidates], n_add)
        hiv.circumcise(new_circs)
        return len(new_circs)

    def step(self):
        sim = self.sim
        ppl = sim.people
        hiv = sim.diseases.hiv

        # Coverage is a circumcision *prevalence* (stock) target, matching
        # cross-sectional survey data: the denominator is ALL eligible males
        # (already-circumcised included), and each step we top up to the target
        # rather than circumcising a fraction of the remaining uncircumcised
        # (which would ratchet coverage toward 100% regardless of the target).
        # Stratified coverage is corrected per (age, sex) stratum so age/sex
        # differentials in the input data are preserved.
        pool = (ppl.male & ppl.alive) & self.check_eligibility()

        n_new = 0
        if self.coverage is not None:
            stratum_targets = compute_stratum_targets(
                self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
                self.ti, pool.uids, sim,
            )
            if stratum_targets is not None:
                for key, target in stratum_targets.items():
                    ab, sex = (key[0], key[1]) if isinstance(key, tuple) else (key, 1)
                    n_new += self._circumcise_to_target(hiv, pool & age_sex_mask(ab, sex, ppl), target)
            else:
                total = compute_coverage_target(
                    self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
                    self.ti, pool.uids, sim,
                )
                if total is not None:
                    n_new += self._circumcise_to_target(hiv, pool, total)

        # Traditional (non-program) circumcision: independent of the coverage-target
        # logic above and not subject to `eligibility`. A one-time probabilistic
        # assessment, made when a pre-debut male first reaches traditional_age;
        # traditional_assessed ensures he's only assessed once even if he doesn't
        # "pass" (so a p=0.1 draw isn't repeated every step he remains eligible).
        if self.pars.traditional_prob.pars.p > 0:
            net = sim.networks.structuredsexual
            newly_eligible = (ppl.male & ppl.alive & ~hiv.circumcised & ~self.traditional_assessed &
                              ~net.over_debut & (ppl.age >= self.pars.traditional_age)).uids
            if len(newly_eligible):
                self.traditional_assessed[newly_eligible] = True
                new_trad = self.pars.traditional_prob.filter(newly_eligible)
                if len(new_trad):
                    hiv.circumcise(new_trad, traditional=True)
                    n_new += len(new_trad)

        self.results['new_circumcisions'][self.ti] = n_new
        self.results['n_circumcised'][self.ti] = count(hiv.circumcised)

        return


class Prep(ss.Intervention):
    """
    Pre-exposure prophylaxis (PrEP).

    Any PrEP product (oral or injectable) is specified as a combination of 
    three parameters: efficacy (``prep_eff``), adherence (``prep_adh``), 
    and duration (``prep_dur``). ``prep_adh`` is a continuous multiplier 
    which reduces the base efficacy (``prep_eff``), ranges from [0-1],
    default value is 1.0.

    The state of the intervention (``on_prep``/``prep_eff``/``prep_source``/`
    `ti_prep_start``/``ti_prep_stop``) lives on :class:`~stisim.diseases.hiv.HIV`.
    This excludes people from enrolling on multiple PrEP products at once.
    
    The ``eligibility`` callable defines the target population. It receives
    the ``Sim`` object and should return a boolean array over all agents.
    HIV-negative and not-already-on-PrEP filters are always applied on top,
    so eligibility only needs to express *who to target*, not the clinical
    preconditions.

    Args:
        prep_eff:         base efficacy (0-1) when fully adherent (default 0.8)
        prep_dur:         course duration before renewal is needed, an
                          ``ss.dur``/Time value (default ``ss.months(3)``)
        prep_adh:         adherence level (0-1), multiplied into prep_eff to get
                          the realized efficacy (default 1.0 = fully adherent)
        coverage:         coverage target; any format accepted by parse_coverage.
                          Defaults to a modest historical ramp-up.
        eligibility:      callable ``(sim) -> BoolArr`` defining the target population;
                          defaults to FSW (``sim.networks.structuredsexual.fsw``)
        smoothness:       interpolation smoothness (0=linear, default)
        format_priority:  when both n_prep and p_prep are non-NaN, prefer this format

    Examples::

        # Default: oral-like PrEP for FSW, ramping up over time
        prep = sti.Prep(coverage={'year': [2020, 2025], 'value': [0, 0.5]})

        # AGYW targeting
        prep = sti.Prep(
            coverage=0.4,
            eligibility=lambda sim: sim.people.female & (sim.people.age < 25),
        )

        # A long-acting product alongside oral PrEP (distinct names required)
        interventions = [
            sti.Prep(prep_eff=0.75, prep_dur=ss.months(3), prep_adh=1.0,
                     coverage=0.3, name='prep_oral'),
                     # source: doi/full/10.1056/NEJMoa1108524
            sti.Prep(prep_eff=0.95, prep_dur=ss.months(12), prep_adh=1.0,
                     coverage=0.1, name='prep_cab_la'),
                     # source: doi.org/10.1016/S0140-6736(22)00538-4
        ]
    """

    @staticmethod
    def _default_eligibility(sim):
        return sim.networks.structuredsexual.fsw

    _source_counter = 0  # class-level counter; each instance gets a unique id

    def __init__(self, pars=None, prep_eff=0.8, prep_dur=None, prep_adh=1.0,
                 coverage=None, eligibility=None, smoothness=0, format_priority='n', **kwargs):
        super().__init__()

        Prep._source_counter += 1
        self._source_id = Prep._source_counter

        self.define_pars(
            prep_eff=prep_eff,
            prep_dur=prep_dur if prep_dur is not None else ss.months(3),
            prep_adh=prep_adh,
        )
        # Pop legacy years= before update_pars rejects it as an unrecognised par
        years = kwargs.pop('years', None)
        if years is not None and coverage is not None:
            coverage = {'year': years, 'value': coverage}
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility if eligibility is not None else self._default_eligibility
        self._smoothness      = smoothness
        self._format_priority = format_priority
        self.coverage         = None
        self.coverage_format  = None
        self.age_bins         = None
        self.sex_keys         = None

        if coverage is None:
            coverage = {'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}
        self._raw_coverage = coverage

        # States
        self.define_states(
            ss.FloatArr('willingness', default=ss.random())
        )

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.coverage, self.coverage_format, self.age_bins, self.sex_keys = parse_coverage(
            self._raw_coverage, valid_names=['n_prep', 'p_prep'], yearvec=self.t.yearvec,
            smoothness=self._smoothness, format_priority=self._format_priority,
        )
        return

    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_prep',  dtype=int, label='New PrEP starts', auto_plot=False),
            ss.Result('n_on_prep', dtype=int, label='Number on PrEP',  auto_plot=False),
        )
        return

    def _add_to_target(self, hiv, pool, target):
        """Top ``pool`` up to ``target`` on THIS instance's product, choosing
        the highest-willingness people not already on any PrEP. Returns the
        number of new starts."""
        on_mine = hiv.on_prep & (hiv.prep_source == self._source_id)
        n_add = int(target) - (pool & on_mine).count()
        if n_add <= 0:
            return 0
        candidates = (pool & ~hiv.on_prep).uids  # excludes ANY product -- mutual exclusivity
        if len(candidates) == 0:
            return 0
        new_uids = top_n_by(candidates, self.willingness[candidates], n_add)
        started = hiv.start_prep(new_uids, eff=self.pars.prep_eff, dur=self.pars.prep_dur,
                                  source_id=self._source_id, adh=self.pars.prep_adh)
        return len(started)

    def step(self):
        sim = self.sim
        ppl = sim.people
        hiv = sim.diseases.hiv

        # Expire courses whose duration has elapsed, for THIS instance's
        # enrollees only -- each Prep instance is responsible for its own.
        expiring = hiv.on_prep & (hiv.prep_source == self._source_id) & (hiv.ti_prep_stop <= self.ti)
        if expiring.any():
            hiv.stop_prep(expiring.uids)

        # Coverage is a PrEP *prevalence* (stock) target, matching how PrEP
        # programs report coverage ("X% of AGYW on PrEP"): the denominator is
        # ALL eligible people (already on this product included), and each step
        # we top up to the target rather than adding a fraction of the
        # remaining pool (which would ratchet toward 100% regardless of the
        # target). Stratified coverage is corrected per (age, sex) stratum.
        pool = ppl.alive & ~hiv.infected & self.eligibility(sim)

        n_new = 0
        if self.coverage is not None:
            stratum_targets = compute_stratum_targets(
                self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
                self.ti, pool.uids, sim,
            )
            if stratum_targets is not None:
                for key, target in stratum_targets.items():
                    # Age-only stratification (no Gender column) means "don't filter by sex" —
                    # the pool's own eligibility already determines who's targeted. VMMC copies
                    # this same pattern but defaults to sex=1 (male), which is wrong here since
                    # Prep's default/typical eligibility (FSW, AGYW) is female.
                    ab, sex = (key[0], key[1]) if isinstance(key, tuple) else (key, None)
                    n_new += self._add_to_target(hiv, pool & age_sex_mask(ab, sex, ppl), target)
            else:
                total = compute_coverage_target(
                    self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
                    self.ti, pool.uids, sim,
                )
                if total is not None:
                    n_new += self._add_to_target(hiv, pool, total)

        self.results['new_prep'][self.ti] = n_new
        self.results['n_on_prep'][self.ti] = (hiv.on_prep & (hiv.prep_source == self._source_id)).count()

        return


