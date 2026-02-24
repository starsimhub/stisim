"""
Tests for hivsim_examples package and hs.demo()
"""
import numpy as np
import hivsim as hs


kw = dict(n_agents=200, dur=5)


def test_simple():
    """Test the simple example via make_sim and hs.demo produce equivalent results."""
    from hivsim_examples.simple.sim import make_sim

    seed = 42
    sim1 = make_sim(rand_seed=seed, **kw)
    sim1.run()

    sim2 = hs.demo('simple', run=False, rand_seed=seed, **kw)
    sim2.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    assert np.allclose(prev1, prev2), 'make_sim and hs.demo should produce identical results'


def test_zimbabwe():
    """Test the Zimbabwe example via make_sim and hs.demo produce equivalent results, and plot."""
    from hivsim_examples.zimbabwe.sim import make_sim

    seed = 42
    sim1 = make_sim(rand_seed=seed, n_agents=200, sim_pars=dict(stop=1995))
    sim1.run()

    sim2 = hs.demo('zimbabwe', run=False, rand_seed=seed, n_agents=200, sim_pars=dict(stop=1995))
    sim2.run()

    prev1 = sim1.results.hiv.prevalence[:]
    prev2 = sim2.results.hiv.prevalence[:]
    assert np.allclose(prev1, prev2), 'make_sim and hs.demo should produce identical results'

    sim1.plot('hiv')


if __name__ == '__main__':
    test_simple()
    test_zimbabwe()
    print('All tests passed!')
