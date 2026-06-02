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

    # ---- Subclass extension points -----------------------------------------
    def _get_pool(self):
        """Eligible agents the kernel operates on. Returns a BoolArr-shaped mask.

        Default: post-debut males. Subclasses override for other pools
        (e.g. heterosexual: return ``self.over_debut`` and add a cross-sex
        constraint in ``_mix_weights_row``).
        """
        return self.over_debut & self.sim.people.male

    def _mix_node_arrays(self):
        """Return a dict of per-node arrays consumed by ``_mix_weights_row``.

        Base implementation provides:
            - ``'nodes'`` (np.ndarray of int uids in the pool, monotonically
              sorted by uid)
            - ``'log1p_deg'`` (1 + log1p of current degree per pool node)

        Subclasses extend the dict with ``'age_vec'``, ``'risk_vec'``, etc.
        as needed and reference them from a custom ``_mix_weights_row``.
        """
        pool_uids = np.asarray(self._get_pool().uids, dtype=np.int64)
        n = pool_uids.size
        if n == 0:
            return {
                'nodes': pool_uids,
                'log1p_deg': np.empty(0, dtype=float),
            }
        # Vectorised degree count: pool_uids is sorted, so searchsorted maps
        # edge endpoints to pool indices in O((|E|+n) log n) without Python
        # loops. Endpoints outside the pool are filtered before bincount.
        endpoints = np.concatenate([
            np.asarray(self.edges.p1, dtype=np.int64),
            np.asarray(self.edges.p2, dtype=np.int64),
        ])
        if endpoints.size:
            idx = np.searchsorted(pool_uids, endpoints)
            in_bounds = idx < n
            hit = np.zeros_like(endpoints, dtype=bool)
            hit[in_bounds] = pool_uids[idx[in_bounds]] == endpoints[in_bounds]
            deg_vec = np.bincount(idx[hit], minlength=n).astype(float)
        else:
            deg_vec = np.zeros(n, dtype=float)
        return {
            'nodes': pool_uids,
            'log1p_deg': 1.0 + np.log1p(deg_vec),
        }

    def _mix_weights_row(self, i, mix_arrays):
        """Vectorised pair weights between node i and all j > i.

        Returns a 1D array of length (n - 1 - i) of non-negative weights.
        Default = preferential-attachment degree multiplier. Subclasses
        may read ``mix_arrays['nodes']`` for per-row uids when needed
        (e.g. to look up an external per-agent state).
        """
        log1p_deg = mix_arrays['log1p_deg']
        return log1p_deg[i] * log1p_deg[i+1:]
