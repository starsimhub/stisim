"""
Placeholder modules for non-specific STIs, including:
    - Genital ulcer disease
    - Vaginal discharge
"""

import starsim as ss
import stisim as sti

__all__ = ['Placeholder', 'GUDPlaceholder', 'VDPlaceholder', 'GUD']


class Placeholder(ss.Disease):
    # A simple placeholder module to use when testing connectors

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
    # Background prevalence of vaginal discharge
    def __init__(self, pars=None, name='vd', **kwargs):
        super().__init__(pars=pars, name=name, **kwargs)
        return


class GUDPlaceholder(Placeholder):
    # Background prevalence of genital ulcerative disease
    def __init__(self, pars=None, name='gud', **kwargs):
        super().__init__(pars=pars, name=name, **kwargs)
        return


class GUD(ss.Infection):

    def __init__(self, pars=None, init_prev_data=None, **kwargs):
        super().__init__()
        self.define_pars(
            dur_inf = ss.lognorm_ex(ss.months(3), ss.months(1)),
            beta=1.0,  # Placeholder
            init_prev=0,  # See make_init_prev_fn
            rel_init_prev=1,
        )
        self.update_pars(pars, **kwargs)

        # Set initial prevalence
        self.init_prev_data = init_prev_data
        if init_prev_data is not None:
            self.pars.init_prev = ss.bernoulli(self.make_init_prev_fn)

        # Add states
        self.define_states(
            ss.FloatArr('ti_recovered'),
        )
        return

    @staticmethod
    def make_init_prev_fn(self, sim, uids):
        return sti.make_init_prev_fn(self, sim, uids)

    def step_state(self):
        """ Update states """
        self.rel_sus[:] = 1
        self.rel_trans[:] = 1

        recovered = (self.infected & (self.ti_recovered <= self.sim.ti)).uids
        self.infected[recovered] = False
        self.susceptible[recovered] = True
        return

    def set_prognoses(self, uids, source_uids=None):
        """
        Set initial prognoses for adults newly infected with syphilis
        """
        ti = self.ti

        self.susceptible[uids] = False
        self.infected[uids] = True
        self.ti_infected[uids] = ti

        # Set future recovery
        dur_inf = self.pars.dur_inf.rvs(uids)
        self.ti_recovered[uids] = ti + dur_inf

        return
