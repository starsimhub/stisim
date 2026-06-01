"""Tests for ``sti.MSMScaleFreeNetwork`` and its FastSampler helper."""
import numpy as np
import pytest
import sciris as sc
import starsim as ss
import stisim as sti
from stisim.networks._fast_sampler import FastSampler


# %% FastSampler unit tests

def test_fast_sampler_uniform_recovery():
    """Uniform weights → uniform samples (within tolerance)."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.ones(4))
    n_draws = 10_000
    counts = np.bincount([sampler.sample_index(rng) for _ in range(n_draws)], minlength=4)
    rel_freq = counts / n_draws
    assert np.all(np.abs(rel_freq - 0.25) < 0.02), f'non-uniform rel_freq: {rel_freq}'


def test_fast_sampler_weighted_recovery():
    """Weight ratio is preserved in the empirical sample."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.array([1.0, 3.0]))
    n_draws = 10_000
    counts = np.bincount([sampler.sample_index(rng) for _ in range(n_draws)], minlength=2)
    rel_freq = counts / n_draws
    # Expected: index 0 -> 0.25, index 1 -> 0.75
    assert abs(rel_freq[0] - 0.25) < 0.02
    assert abs(rel_freq[1] - 0.75) < 0.02


def test_fast_sampler_zero_mass_raises():
    """Empty or all-zero input raises on sample."""
    rng = np.random.default_rng(0)
    for w in [np.array([]), np.zeros(5)]:
        sampler = FastSampler(w)
        with pytest.raises(RuntimeError, match='no mass'):
            sampler.sample_index(rng)


def test_fast_sampler_handles_negative_weights():
    """Negative weights are clipped to zero (not raise)."""
    rng = np.random.default_rng(0)
    sampler = FastSampler(np.array([-1.0, 2.0, -3.0, 4.0]))
    # Only indices 1 and 3 should ever be drawn
    drawn = {sampler.sample_index(rng) for _ in range(200)}
    assert drawn <= {1, 3}, f'drew clipped-out indices: {drawn}'
