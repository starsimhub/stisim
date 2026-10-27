"""
Genital ulcer disease (GUD), and placeholder modules for non-specific STIs (genital ulcer disease and vaginal discharge) that keep a target prevalence
"""

import numpy as np
import starsim as ss
import stisim as sti

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


class GUD(ss.Infection):
    """
    Genital ulcer disease (GUD): a simple SIS infection, sexually transmitted, with a lognormal duration

    Typically used alongside a syphilis connector that drives GUD prevalence from active syphilis cases.

    Args:
        pars (dict): parameter values, as for the keywords
        name (str): the name of the module (default ``'gud'``)
        kwargs (dict): parameter values (``beta``, ``dur_inf``), ``init_prev``, or ``label``
    """
    def __init__(self, pars=None, name='gud', **kwargs):
        super().__init__(name=name, label=kwargs.pop('label', None))
        self.define_pars(
            beta = 1.0, # Probability of transmission per act on the sexual networks (a placeholder)
            dur_inf = ss.lognorm_ex(ss.months(3), ss.months(1)), # Duration of infection
            init_prev = 0,
        )
        self.update_pars(pars, **kwargs)
        self.define_states(ss.TimeArr('ti_clear', label='Time of clearance'))
        return

    def validate_beta(self):
        """ The same beta on each sexual network, both ways, and none on the others """
        beta = self.pars.beta
        self.pars.beta = sti.sexual_beta(beta, beta, beta)
        try:
            return super().validate_beta()
        finally:
            self.pars.beta = beta

    def set_prognoses(self, uids, sources=None):
        super().set_prognoses(uids, sources)
        self.infected[uids] = True
        self.ti_infected[uids] = self.ti if self.sim.initialized else 0 # The infections seeded at the start are of its first point (see SEIS.set_prognoses())
        start = self.ti + 1 if self.sim.initialized else self.ti # A course runs from when the infection happens: at the end of the step, in transmission (as a declared disease dates it), or as ss.Infection.init_post() dates those seeded at the start
        self.ti_clear[uids] = start + self.pars.dur_inf.rvs(uids)
        return

    def step_state(self):
        """ Clearance, back to susceptible """
        uids = (self.infected & self.ti_clear.due()).uids
        self.infected[uids] = False
        self.susceptible[uids] = True
        self.ti_clear[uids] = np.nan
        return
