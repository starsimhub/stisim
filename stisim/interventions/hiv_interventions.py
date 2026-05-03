"""
Define HIV interventions for STIsim
By default, these all have units of a year and timesteps of 1/12
"""

import starsim as ss
import numpy as np
from stisim.interventions.base_interventions import STITest
from stisim.interventions.utils import (
    parse_coverage, compute_coverage_target,
)
from stisim.utils import count


# %% HIV classes
__all__ = ["HIVDx", "HIVTest", "InfantHIVTest", "ART", "VMMC", "Prep"]


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
        3. ART checks for newly diagnosed agents (ti_diagnosed == current ti)
        4. Newly diagnosed agents initiate ART with probability art_initiation
        5. If coverage data is provided, ART corrects to match targets

    Args:
        test_prob_data: annual testing probability (if dt_scale=True, the default).
            A value of 0.1 means ~10% of eligible agents tested per year. To
            specify a per-timestep probability instead, set dt_scale=False.
        eligibility (func): who can be tested. Default: undiagnosed agents.
        start (float): calendar year when testing begins
        dt_scale (bool): if True (default), test_prob_data is an annual probability.
            Set to False for per-timestep probability.

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
            eligibility=lambda sim: sim.demographics.pregnancy.tri1_uids[
                ~sim.diseases.hiv.diagnosed[sim.demographics.pregnancy.tri1_uids]
            ],
        )
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None,
                 eligibility=None, name=None, label=None, newborn_test=None, **kwargs):
        if product is None: product = HIVDx(name=f'HIVDx_{name}')
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years,
                         start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed
        self.newborn_test = newborn_test

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = self.ti

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
            ti_births = pregnancy.ti_delivery[mn.p1[pos_mask]]
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


class ART(ss.Intervention):
    """
    Antiretroviral therapy intervention.

    Requires HIVTest (or equivalent) to diagnose agents first — ART only
    initiates agents who have hiv.diagnosed=True. A warning is raised if no
    HIVTest is found in the sim.

    Processing flow each timestep:
        1. Agents scheduled to stop ART are removed
        2. Newly diagnosed agents (ti_diagnosed == this timestep) are filtered
           by art_initiation probability
        3. If coverage is specified: agents are added/removed to match the target
           number, prioritized by CD4 count and care-seeking propensity
        4. If no coverage is specified: all newly diagnosed who pass art_initiation
           go directly on ART (no capacity constraint)
        5. Mothers on ART reduce infant susceptibility (prenatal via MaternalNet,
           postnatal via BreastfeedingNet) by pmtct_efficacy

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
        smoothness:       interpolation smoothness (0=linear, default)
        format_priority:  when both n_art and p_art are non-NaN, prefer this format
                          ('n' or 'p', default 'n')

    Examples::

        # Simple: 80% of infected on ART
        art = sti.ART(coverage=0.8)

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

    def __init__(self, pars=None, coverage=None, smoothness=0, format_priority='n', **kwargs):
        super().__init__()

        self.define_pars(
            art_initiation=ss.bernoulli(p=0.9),
            pmtct_efficacy=0.96,  # How much maternal ART reduces infant susceptibility
                                   # to HIV via MaternalNet (prenatal) and BreastfeedingNet
                                   # (postnatal). Conceptually like infant PrEP.
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._smoothness      = smoothness
        self._format_priority = format_priority
        self.coverage         = None  # Set in init_pre
        self.coverage_format  = None  # 'n', 'p', or per-timestep array
        self.age_bins         = None  # For stratified coverage
        self.sex_keys         = None  # For stratified coverage
        return

    def init_pre(self, sim):
        super().init_pre(sim)

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
        self.initialized = True
        return

    def _get_n_to_treat(self, eligible_uids):
        """Get the target number of people on ART this timestep."""
        return compute_coverage_target(
            self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
            self.ti, eligible_uids, self.sim,
        )

    def step(self):
        """
        Apply ART at each timestep: stop ART for those scheduled, initiate for newly diagnosed,
        and correct overall coverage to match targets.
        """
        sim = self.sim
        hiv = sim.diseases.hiv
        inf_uids = hiv.infected.uids

        # Determine treatment target (None = no capacity constraint)
        n_to_treat = self._get_n_to_treat(inf_uids)

        # Check who is stopping ART
        if hiv.on_art.any():
            stopping = hiv.on_art & (hiv.ti_stop_art <= self.ti)
            if stopping.any():
                try:
                    hiv.stop_art(stopping.uids)
                except:
                    errormsg = f'Error stopping ART for {stopping.uids}'
                    raise ValueError(errormsg)

        # Initiate ART for newly diagnosed
        diagnosed = hiv.ti_diagnosed == self.ti
        if len(diagnosed.uids):
            dx_to_treat = self.pars.art_initiation.filter(diagnosed.uids)

            if n_to_treat is None:
                # No coverage target — treat all who initiate
                hiv.start_art(dx_to_treat)
            else:
                # Coverage target — only treat if spots available
                on_art = hiv.on_art
                n_available_spots = n_to_treat - len(on_art.uids)
                if n_available_spots > 0:
                    self.prioritize_art(sim, n=n_available_spots, awaiting_art_uids=dx_to_treat)

        # Correct coverage to match target (only if target is set)
        if n_to_treat is not None:
            self.art_coverage_correction(sim, target_coverage=n_to_treat)

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
                art_pregnant = (hiv.on_art & preg.pregnant).uids
                if len(art_pregnant):
                    unborn = sim.networks.maternalnet.find_contacts(art_pregnant)
                    hiv.rel_sus[ss.uids(unborn)] *= 1 - pmtct_eff

            # Postnatal: protect breastfed infants of breastfeeding mothers on ART
            if hasattr(sim.networks, 'breastfeedingnet'):
                art_breastfeeding = (hiv.on_art & preg.breastfeeding).uids
                if len(art_breastfeeding):
                    breastfed = sim.networks.breastfeedingnet.find_contacts(art_breastfeeding)
                    hiv.rel_sus[ss.uids(breastfed)] *= 1 - pmtct_eff

        return

    def prioritize_art(self, sim, n=None, awaiting_art_uids=None):
        """ Prioritize ART to n agents among those awaiting treatment """
        hiv = sim.diseases.hiv
        if awaiting_art_uids is None:
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids

        # Enough spots for everyone
        if n > len(awaiting_art_uids):
            start_uids = awaiting_art_uids

        # Not enough spots — prioritize by CD4 and care seeking
        else:
            cd4_counts   = hiv.cd4[awaiting_art_uids]
            care_seeking = hiv.care_seeking[awaiting_art_uids]
            weights = cd4_counts * (1 / care_seeking)
            choices = np.argsort(weights)[:n]
            start_uids = awaiting_art_uids[choices]

        hiv.start_art(start_uids)

        return

    def art_coverage_correction(self, sim, target_coverage=None):
        """ Adjust ART coverage to match data """
        hiv = sim.diseases.hiv
        on_art = hiv.on_art

        # Too many on treatment → remove
        if len(on_art.uids) > target_coverage:
            n_to_stop    = int(len(on_art.uids) - target_coverage)
            on_art_uids  = on_art.uids
            cd4_counts   = hiv.cd4[on_art_uids]
            care_seeking = hiv.care_seeking[on_art_uids]
            weights  = cd4_counts / care_seeking
            choices  = np.argsort(-weights)[:n_to_stop]
            stop_uids = on_art_uids[choices]
            hiv.ti_stop_art[stop_uids] = self.ti
            hiv.stop_art(stop_uids)

        # Not enough on treatment → add
        elif len(on_art.uids) < target_coverage:
            n_to_add = target_coverage - len(on_art.uids)
            awaiting_art_uids = (hiv.diagnosed & ~hiv.on_art).uids
            self.prioritize_art(sim, n=n_to_add, awaiting_art_uids=awaiting_art_uids)

        return


class VMMC(ss.Intervention):
    """
    Voluntary medical male circumcision.

    Reduces male susceptibility to HIV acquisition by eff_circ (default 60%).
    Unlike ART, VMMC does not require diagnosis — it circumcises males up to a
    coverage target, prioritized by willingness (a random per-agent score).

    If no coverage is specified, VMMC does nothing. Coverage must be provided
    explicitly via the coverage parameter.

    Coverage is parsed by :func:`parse_coverage` (same formats as ART, using
    ``n_vmmc``/``p_vmmc`` column names). Age-only stratification is supported
    (no Gender column required, since VMMC is males-only).

    Args:
        coverage:         coverage target (default None; VMMC does nothing without data).
                          See :func:`parse_coverage` for supported formats.
        eff_circ:         efficacy (default 0.6 = 60% reduction in HIV acquisition)
        eligibility:      optional function to restrict who is eligible (default: all males)
        smoothness:       interpolation smoothness (0=linear, default)
        format_priority:  when both n_vmmc and p_vmmc are non-NaN, prefer this format

    Examples::

        vmmc = sti.VMMC(coverage=0.3)
        vmmc = sti.VMMC(coverage={'year': [2010, 2025], 'value': [0, 0.4]})
    """

    def __init__(self, pars=None, coverage=None, eligibility=None, smoothness=0, format_priority='n', **kwargs):
        super().__init__(eligibility=eligibility)

        self.define_pars(
            eff_circ=0.6,
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
        self.willingness     = ss.FloatArr('willingness', default=ss.random())
        self.circumcised     = ss.BoolArr('circumcised', default=False)
        self.ti_circumcised  = ss.FloatArr('ti_circumcised')

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

    def _get_n_to_circ(self, eligible_uids):
        """Get the target number of circumcisions this timestep."""
        return compute_coverage_target(
            self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
            self.ti, eligible_uids, self.sim,
        )

    def step(self):
        sim = self.sim
        
        # Get eligible people by combining user-specified eligibility function with male & uncircumcised
        eligible_uids = self.check_eligibility()
        eligible_uids = ((sim.people.male & ~self.circumcised) & eligible_uids).uids

        n_to_circ = self._get_n_to_circ(eligible_uids)

        if n_to_circ is not None and n_to_circ > 0:
            weights = self.willingness[eligible_uids]
            choices = np.argsort(-weights)[:n_to_circ]
            new_circs = eligible_uids[choices]

            self.circumcised[new_circs] = True
            self.ti_circumcised[new_circs] = self.ti

        self.results['new_circumcisions'][self.ti] = n_to_circ or 0
        self.results['n_circumcised'][self.ti] = count(self.circumcised)

        # Reduce rel_sus
        sim.diseases.hiv.rel_sus[self.circumcised] *= 1 - self.pars.eff_circ

        return


class Prep(ss.Intervention):
    """
    Pre-exposure prophylaxis (PrEP).

    Reduces HIV susceptibility by ``eff_prep`` (default 80%) among eligible agents.
    Coverage ramps up over time via ``parse_coverage`` (same flexible inputs as
    ART/VMMC). Uses a per-agent probability model (coverage = probability of being
    on PrEP each step) rather than a target-count model like ART/VMMC.

    The ``eligibility`` callable defines the target population. It receives the
    ``Sim`` object and should return a boolean array over all agents. HIV-negative
    and not-already-on-PrEP filters are always applied on top, so eligibility only
    needs to express *who to target*, not the clinical preconditions.

    Args:
        coverage:    coverage level(s); any format accepted by parse_coverage,
                     or legacy (years, coverage) list pairs via pars
        eff_prep:    efficacy (default 0.8 = 80% reduction in acquisition)
        smoothness:  interpolation smoothness (0=linear, default)
        eligibility: callable ``(sim) -> BoolArr`` defining the target population;
                     defaults to FSW (``sim.networks.structuredsexual.fsw``)

    Examples::

        # Default: FSW at time-varying coverage
        prep = sti.Prep(coverage={'year': [2020, 2025], 'value': [0, 0.5]})

        # AGYW targeting
        prep = sti.Prep(
            coverage=0.4,
            eligibility=lambda sim: sim.people.female & (sim.people.age < 25),
        )
    """

    @staticmethod
    def _default_eligibility(sim):
        return sim.networks.structuredsexual.fsw

    def __init__(self, pars=None, coverage=None, eligibility=None, smoothness=0, **kwargs):
        super().__init__()
        self.define_pars(
            coverage_dist=ss.bernoulli(p=0),
            eff_prep=0.8,
        )
        # Pop legacy years= before update_pars rejects it as an unrecognised par
        years = kwargs.pop('years', None)
        if years is not None and coverage is not None:
            coverage = {'year': years, 'value': coverage}
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility if eligibility is not None else self._default_eligibility
        self._smoothness = smoothness
        self._coverage_arr = None  # Set in init_pre

        if coverage is None:
            coverage = {'year': [2004, 2005, 2015, 2025], 'value': [0, 0.01, 0.5, 0.8]}
        self._raw_coverage = coverage

        self.define_states(
            ss.BoolArr('on_prep', label='On PrEP'),
        )
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self._coverage_arr, _, _, _ = parse_coverage(
            self._raw_coverage, valid_names=['p_prep'], yearvec=self.t.yearvec,
            smoothness=self._smoothness,
        )
        return

    def step(self):
        sim = self.sim
        cov_val = self._coverage_arr[self.ti]
        if cov_val > 0:
            self.pars.coverage_dist.set(p=cov_val)
            eligible = self.eligibility(sim) & ~sim.diseases.hiv.infected & ~self.on_prep
            new_on_prep = self.pars.coverage_dist.filter(eligible)
            sim.diseases.hiv.rel_sus[new_on_prep] *= 1 - self.pars.eff_prep
        return


