"""
Minimal HIV simulation example using hivsim defaults.
"""
import hivsim


def make_sim(**kwargs):
    """Create a simple HIV simulation with default parameters."""
    kwargs.setdefault('n_agents', 2000)
    kwargs.setdefault('dur', 20)
    return hivsim.Sim(**kwargs)
