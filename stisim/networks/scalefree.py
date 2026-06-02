"""Scale-free / preferential-attachment sexual network.

Implements the Whittles-2019 S2 pair-formation kernel as ``MSMScaleFreeNetwork``,
restricted to post-debut male agents. The kernel is sex-agnostic; subclasses
can override ``_get_pool()`` for other configurations.

Ported from ``starsim_x/BespokeNet`` per the cross-walk (May 2026), restricted
to the reusable L2 primitives.

References:
    Whittles LK et al. (2019) Disease control, demography and the role of
    pre-exposure prophylaxis: a modelling study among MSM in England.
"""
import numpy as np
import starsim as ss

from .base import BaseNetwork
from ._fast_sampler import FastSampler


__all__ = ['MSMScaleFreeNetwork']


class MSMScaleFreeNetwork(BaseNetwork):
    """Preferential-attachment MSM sexual network with continuous-time
    formation and Markovian deletion.

    Eligible agents are post-debut males. Edges form via a rate kernel
    proportional to ``g(λ_i) g(λ_j) × pair_weight``, where the default
    pair weight is a degree-based multiplier (rich-get-richer). Edges are
    deleted uniformly at random subject to a hard duration cap.

    Subclasses can plug in age- or risk-weighted formation by overriding
    ``_mix_node_arrays`` and ``_mix_weights_row``.

    Args:
        target_mean_degree (float): target mean number of concurrent
            partners per agent in the pool (default 2.0).
        target_mean_dur (ss.dur or int): target mean edge duration. Accepts
            an ``ss.dur`` (converted to integer steps at ``init_pre`` via
            ``self.t.dt``) or a raw integer step count. Default ``ss.years(2)``.
        max_edge_dur (ss.dur or int): hard cap on edge persistence. Default
            ``ss.years(10)``.
        phi (float): Whittles-2019 turnover parameter. Default 1.0. Sets
            the q0 initial-network density via ``q0 = 1/(1+phi)``.
        name (str): network name (default: auto-assigned).
        **kwargs: forwarded to ``update_pars``.
    """

    def __init__(self, pars=None, name=None, **kwargs):
        super().__init__(name=name)
        self.define_pars(
            target_mean_degree=2.0,
            target_mean_dur=ss.years(2),
            max_edge_dur=ss.years(10),
            phi=1.0,
        )
        self.update_pars(pars, **kwargs)
        # Internal step-domain copies of the user-facing dur parameters; set
        # by init_pre once self.t.dt is known.
        self._target_mean_dur_steps = None
        self._max_edge_dur_steps = None
        # Kernel artefacts cached between rebuilds; populated by _build_kernel.
        self._kernel_pairs_i = None
        self._kernel_pairs_j = None
        self._kernel_add_sampler = None
        self._kernel_hat_Ra = None
        self._kernel_rebuild_count = 0
        # Hardcoded BespokeNet defaults — promote to user pars only when needed.
        self._k1 = 1.0
        self._g_normalization_guard = 0.25
        self._rebuild_every = 10
        self._age_a = 0.0
        self._partnering_assortivity = 1.0
        return

    def init_pre(self, sim):
        """Convert ``ss.dur`` parameters to integer step counts.

        Accepts an int directly (treated as a raw step count) so users can
        bypass ``ss.dur`` if they prefer.
        """
        super().init_pre(sim)
        self._target_mean_dur_steps = self._as_steps(self.pars.target_mean_dur)
        self._max_edge_dur_steps = self._as_steps(self.pars.max_edge_dur)
        if self._max_edge_dur_steps < self._target_mean_dur_steps:
            raise ValueError(
                f'max_edge_dur ({self._max_edge_dur_steps} steps) must be >= '
                f'target_mean_dur ({self._target_mean_dur_steps} steps)'
            )
        return

    def _as_steps(self, dur):
        """Convert an ``ss.dur`` to integer step counts using ``self.t.dt``.

        Raw ints are returned unchanged.
        """
        if isinstance(dur, (int, np.integer)):
            return int(dur)
        return int(round(float(dur) / float(self.t.dt)))
