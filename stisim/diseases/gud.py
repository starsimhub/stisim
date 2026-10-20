"""
Placeholder modules for non-specific STIs, including:
    - Genital ulcer disease
    - Vaginal discharge
"""

import starsim as ss
import stisim as sti
from stisim.diseases.sti import BaseSTI

__all__ = ['Placeholder', 'GUDPlaceholder', 'VDPlaceholder', 'GUD']


class Placeholder(ss.Disease):
    """
    Simple placeholder disease module for testing connectors.

    Maintains a target prevalence by probabilistically adding or removing
    symptomatic agents each timestep. Useful for testing disease interaction
    connectors without a full transmission model.

    Args:
        pars (dict): Override default parameters (prevalence, care_seeking).
        name (str): Module name.
        **kwargs: Additional parameters passed to ``update_pars``.
    """

    def __init__(self, pars=None, name=None, **kwargs):
        super().__init__(name=name)

        self.define_pars(
            prevalence=0.1,  # Target prevalance. If None, no automatic infections will be applied
            care_seeking=ss.bernoulli(p=0.5),
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('symptomatic'),  # Symptomatic
            ss.FloatArr('ti_symptomatic'),  # Time of active symptoms
            ss.BoolState('seeking_care'),  # Care seeking
            ss.FloatArr('ti_seeks_care'),  # Time of active symptoms
        )
        self._prev_dist = ss.bernoulli(p=0)

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if not isinstance(self.pars.prevalence, sti.TimeSeries):
            ts = sti.TimeSeries(assumption=self.pars.prevalence)
        else:
            ts = self.pars.prevalence
        self._target_prevalence = ts.interpolate(sim.timevec)

    def set_prognoses(self, target_uids, source_uids=None):
        self.symptomatic[target_uids] = True

    def step(self):
        """
        When using a connector to the syphilis module, this is not needed. The connector should update the syphilis-positive state.
        """

        if self.pars.prevalence is None:
            return

        sim = self.sim

        # Get current prevalence
        n_symptomatic = self.symptomatic.count()
        prev = n_symptomatic/len(sim.people)
        target = self._target_prevalence[sim.ti]
        change = target-prev

        if change > 0:
            # Add a proportion of people that are not infected
            uids = self.symptomatic.false()
            self._prev_dist.set(p=change/(len(uids)/len(sim.people)))
            self.symptomatic[self._prev_dist.filter(uids)] = True
        elif change < 0:
            uids = self.symptomatic.true()
            self._prev_dist.set(p=-change/(len(uids)/len(sim.people)))
            self.symptomatic[self._prev_dist.filter(uids)] = False


class VDPlaceholder(Placeholder):
    """Placeholder for background prevalence of vaginal discharge."""
    def __init__(self, pars=None, name='vd', **kwargs):
        super().__init__(pars=pars, name=name, **kwargs)
        return


class GUDPlaceholder(Placeholder):
    """Placeholder for background prevalence of genital ulcerative disease."""
    def __init__(self, pars=None, name='gud', **kwargs):
        super().__init__(pars=pars, name=name, **kwargs)
        return


class GUD(BaseSTI):
    """
    Genital ulcer disease (GUD), as an SIS infection

    States: S (susceptible) and I (an ulcer), which clears after `dur_inf`. Transmission is
    per act over the sexual network (`beta_m2f`, `rel_beta_f2m`), as for the other STIs; the
    defaults are placeholders, and with `init_prev=0` nobody has GUD unless it is seeded.
    Usually used with a syphilis connector (`sti.gud_syph`).

    Args:
        init_prev (float): initial prevalence (default 0)
        kwargs (dict): passed to `sti.BaseSTI` (e.g. `dur_inf=`, `beta_m2f=`)

    Examples:
        ```python
        gud = sti.GUD(init_prev=0.05, beta_m2f=0.2)
        sim = sti.Sim(diseases=[gud, 'syph'], connectors=sti.gud_syph(gud, 'syph'))
        ```
    """
    S, I = ss.states('S I')
    dur_inf = ss.par(ss.lognorm_ex(ss.months(3), ss.months(1)), doc='Duration of an ulcer')
    infection = S >> I @ BaseSTI.beta*I
    recovery = ss.Transition(I, S, dur=dur_inf)
    n_infected = I
    prevalence = I/(S + I)
    incidence = infection/S
    defaults = dict(beta_m2f=1.0, rel_beta_f2m=1.0, init_prev=0) # Starsim 3: beta=1 on every network (a placeholder)

    def __init__(self, name='gud', **kwargs):
        super().__init__(name=name, **kwargs)
        return

    def flat_results(self):
        return super().flat_results() + [('prevalence', 'prevalence', True, 'Prevalence')]

    @property
    def susceptible(self): return self.S
    @property
    def infected(self):    return self.I
