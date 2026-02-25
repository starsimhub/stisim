"""
Minimal HIV simulation example using hivsim defaults.
"""
import hivsim as hs


def make_sim(**kwargs):
    """Create a simple HIV simulation with default parameters."""
    kwargs.setdefault('n_agents', 2000)
    kwargs.setdefault('dur', 20)
    return hs.Sim(**kwargs)
