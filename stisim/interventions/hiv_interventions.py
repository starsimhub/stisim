"""
Define HIV interventions for STIsim
By default, these all have units of a year and timesteps of 1/12
"""

import starsim as ss
import numpy as np
from stisim.interventions.base_interventions import STITest
from stisim.interventions.utils import (
    parse_coverage, _handle_deprecated_coverage, compute_coverage_target,
)
from stisim.utils import count


# %% HIV classes
__all__ = ["HIVDx", "HIVTest", "ART", "VMMC", "Prep"]


class HIVDx(ss.Product):
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
    """
    def __init__(self, product=None, pars=None, test_prob_data=None, years=None, start=None, eligibility=None, name=None, label=None, **kwargs):
        if product is None: product = HIVDx(name=f'HIVDx_{name}')
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years, start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = self.ti
        return outcomes


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
        5. Mothers on ART protect unborn infants (rel_sus=0) via MaternalNet

    Coverage can be specified in several formats:
        - None: no coverage target; treat all who initiate (default)
        - Scalar (e.g. 0.8): constant proportion of infected on ART
        - Dict: {'year': [2000, 2020], 'value': [0, 0.9]} — linearly interpolated
        - DataFrame: index=years, column 'n_art' (absolute numbers) or 'p_art'
          (proportion of infected). Values are linearly interpolated to the sim
          yearvec.
        - Stratified DataFrame: columns Year, Gender/Sex, AgeBin (format
          '[lo,hi)'), plus a numeric value column. Values are interpolated per
          stratum.

    Intervention ordering: HIVTest must appear before ART in the interventions
    list so that agents diagnosed this timestep can initiate ART in the same step.

    Args:
        coverage:         coverage target in any format above (default None)
        art_initiation:   probability a newly diagnosed person initiates ART
                          (default: ss.bernoulli(p=0.9)). Set to 1 to treat all
                          diagnosed.

    Example::

        # Simple: 80% of infected on ART
        art = sti.ART(coverage=0.8)

        # Time-varying coverage
        art = sti.ART(coverage={'year': [2000, 2010, 2025], 'value': [0, 0.5, 0.9]})

        # No coverage target — 90% of newly diagnosed initiate (default art_initiation)
        art = sti.ART()

        # Treat ALL diagnosed with no coverage constraint
        art = sti.ART(art_initiation=1)

        # From a CSV file
        art = sti.ART(coverage=pd.read_csv('art_coverage.csv').set_index('year'))
    """

    def __init__(self, pars=None, coverage=None, coverage_data=None, smoothness=0, **kwargs):
        super().__init__()

        # Handle deprecated kwargs
        coverage, future_coverage = _handle_deprecated_coverage(coverage, coverage_data, kwargs)

        self.define_pars(
            art_initiation=ss.bernoulli(p=0.9),
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._future_coverage = future_coverage  # Legacy compat, removed once deprecated
        self._smoothness      = smoothness
        self.coverage         = None  # Set in init_pre
        self.coverage_format  = None  # 'n' or 'p'
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
            smoothness=self._smoothness,
        )
        self.initialized = True
        return

    def _get_n_to_treat(self, eligible_uids):
        """Get the target number of people on ART this timestep."""
        return compute_coverage_target(
            self.coverage, self.coverage_format, self.age_bins, self.sex_keys,
            self.ti, eligible_uids, self.sim,
            future_coverage=self._future_coverage,
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

        # Adjust rel_sus for protected unborn agents (only if pregnancy is modeled)
        if hasattr(sim.people, 'pregnancy') and hasattr(sim.networks, 'maternalnet'):
            if hiv.on_art[sim.people.pregnancy.pregnant].any():
                mother_uids = (hiv.on_art & sim.people.pregnancy.pregnant).uids
                infants = sim.networks.maternalnet.find_contacts(mother_uids)
                hiv.rel_sus[ss.uids(infants)] = 0

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

    Coverage formats (same as ART):
        - Scalar: constant proportion of males (e.g. 0.3)
        - Dict: {'year': [...], 'value': [...]} — linearly interpolated
        - DataFrame: index=years, column 'n_vmmc' or 'p_vmmc'
        - Stratified DataFrame: Year/Gender/AgeBin columns

    Args:
        coverage:    coverage target in any format above (default None; VMMC does
                     nothing without coverage data)
        eff_circ:    efficacy (default 0.6 = 60% reduction in HIV acquisition)
        eligibility: optional function to restrict who is eligible (default: all males)

    Example::

        vmmc = sti.VMMC(coverage=0.3)
        vmmc = sti.VMMC(coverage={'year': [2010, 2025], 'value': [0, 0.4]})
    """

    def __init__(self, pars=None, coverage=None, coverage_data=None, eligibility=None, smoothness=0, **kwargs):
        super().__init__(eligibility=eligibility)

        # Handle deprecated kwargs
        coverage, future_coverage = _handle_deprecated_coverage(coverage, coverage_data, kwargs)

        self.define_pars(
            eff_circ=0.6,
        )
        self.update_pars(pars, **kwargs)

        self._raw_coverage    = coverage
        self._future_coverage = future_coverage
        self._smoothness      = smoothness
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
            smoothness=self._smoothness,
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
            future_coverage=self._future_coverage,
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

    By default targets HIV-negative FSWs who are not already on PrEP.
    Reduces HIV susceptibility by eff_prep (default 80%). Coverage ramps up
    over time via ``parse_coverage`` (same flexible inputs as ART/VMMC).
    Use the eligibility parameter to target a different population.

    Note: PrEP uses a per-agent probability model (coverage = probability of
    being on PrEP) rather than a target-count model like ART/VMMC.

    Args:
        coverage:    coverage data in any format accepted by parse_coverage;
                     also accepts legacy (years, coverage) list pairs via pars
        eff_prep:    efficacy (default 0.8 = 80% reduction in acquisition)
        smoothness:  interpolation smoothness (0=linear, default)
        eligibility: function to override default FSW targeting

    Examples::

        prep = sti.Prep(coverage={'year': [2020, 2025], 'value': [0, 0.5]})
        prep = sti.Prep(coverage=0.3)
    """
    def __init__(self, pars=None, coverage=None, eligibility=None, smoothness=0, **kwargs):
        super().__init__()
        self.define_pars(
            coverage_dist=ss.bernoulli(p=0),
            eff_prep=0.8,
        )
        self.update_pars(pars, **kwargs)
        self.eligibility = eligibility
        self._smoothness = smoothness
        self._coverage_arr = None  # Set in init_pre

        # Support legacy (years, coverage) pars by converting to dict format
        if coverage is None and 'years' in self.pars and 'coverage' in self.pars:
            coverage = {'year': self.pars.years, 'value': self.pars.coverage}
        elif coverage is None:
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
            el_fsw = sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.infected & ~self.on_prep
            fsw_on_prep = self.pars.coverage_dist.filter(el_fsw)
            sim.diseases.hiv.rel_sus[fsw_on_prep] *= 1 - self.pars.eff_prep

        return


