"""Inverse-CDF discrete sampler used by ``scalefree.py`` for the S2 kernel.

Build once in O(m); each sample is O(log m). Internal helper — not exported
from ``stisim.networks``. Cross-walk flags this as a starsim L1 upstream
candidate; promote to ``ss.FastSampler`` in a follow-up PR.

Ported from ``starsim_x/networks.py:1315-1354`` (Aug 2025).
"""
import numpy as np


class FastSampler:
    """O(log m) sampling from a discrete weighted distribution.

    Args:
        weights: 1D array of non-negative weights. Need not be normalised.
            Empty or all-zero weights produce a degenerate sampler that
            raises ``RuntimeError`` on sample.
    """

    def __init__(self, weights):
        w = np.asarray(weights, dtype=float).ravel()
        if w.ndim != 1:
            raise ValueError('weights must be 1D')
        if w.size == 0:
            self._cdf = np.array([1.0])
            self._n = 0
            return
        w = np.maximum(w, 0.0)
        total = float(w.sum())
        if total <= 0.0:
            self._cdf = np.array([1.0])
            self._n = 0
            return
        cs = np.cumsum(w)
        self._cdf = cs / total
        self._cdf[-1] = 1.0  # guard against floating-point drift
        self._n = w.size

    def sample_index(self, rng):
        """Return a single index drawn from the weight distribution.

        Args:
            rng: a ``numpy.random.Generator`` instance.

        Raises:
            RuntimeError: if the sampler holds no weight (empty or
                all-zero input).
        """
        if self._n == 0:
            raise RuntimeError('FastSampler: no mass to sample from')
        u = float(rng.random())
        k = int(np.searchsorted(self._cdf, u, side='left'))
        if k >= self._n:
            k = self._n - 1
        return k
