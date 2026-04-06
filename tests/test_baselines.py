"""
Test that the current version of STIsim exactly matches
the baseline results, and benchmark performance.

Uses the Zimbabwe HIV example as the reference simulation.
"""
import sciris as sc
import starsim as ss
import hivsim

baseline_filename  = sc.thisdir(__file__, 'baseline.yaml')
benchmark_filename = sc.thisdir(__file__, 'benchmark.yaml')
sc.options(interactive=False) # Assume not running interactively

# Default parameters for testing (smaller/shorter for speed)
default_kw = dict(n_agents=2000, stop=2010, rand_seed=2, verbose=0)


def make_sim(run=False, **kwargs):
    """
    Create a Zimbabwe HIV simulation for baseline/benchmark testing.
    """
    kw = sc.mergedicts(default_kw, kwargs)
    sim = hivsim.demo('zimbabwe', run=False, **kw)

    if run:
        sim.run()
        sim.plot()

    return sim


def save_baseline():
    """
    Refresh the baseline results. This function is not called during standard
    testing, but instead is called directly when the baseline needs updating.
    """
    sc.heading('Updating baseline values...')

    sim = make_sim()
    sim.run()

    json = sim.to_json(keys='summary')['summary']
    sc.saveyaml(baseline_filename, json)

    print('Done.')
    return


@sc.timer()
def test_baseline():
    """ Compare the current default sim against the saved baseline """

    # Load existing baseline
    old = sc.loadyaml(baseline_filename)

    # Calculate new baseline
    new = make_sim()
    new.run()

    # Compute the comparison
    ss.diff_sims(old, new, die=True)

    return new


@sc.timer()
def test_benchmark(do_save=False, repeats=1, verbose=True):
    """ Compare benchmark performance """

    if verbose: print('Running benchmark...')
    try:
        previous = sc.loadyaml(benchmark_filename)
    except FileNotFoundError:
        previous = None

    t_inits = []
    t_runs  = []
    ref = 270 # Reference benchmark for sc.benchmark(which='numpy') on an Intel i7-12700H (for scaling performance)

    # Test CPU performance before the run
    r1 = sc.benchmark(which='numpy')

    # Do the actual benchmarking
    for r in range(repeats):

        print(f'Repeat {r}')

        # Time initialization
        t0 = sc.tic()
        sim = make_sim()
        sim.init()
        t_init = sc.toc(t0, output=True)

        # Time running
        t0 = sc.tic()
        sim.run()
        t_run = sc.toc(t0, output=True)

        # Store results
        t_inits.append(t_init)
        t_runs.append(t_run)

    # Test CPU performance after the run
    r2 = sc.benchmark(which='numpy')
    ratio = (r1+r2)/2/ref
    t_init = ratio*min(t_inits)
    t_run  = ratio*min(t_runs)

    # Construct json
    n_decimals = 3
    json = {'time': {
                'initialize': round(t_init, n_decimals),
                'run':        round(t_run,  n_decimals),
                },
            'parameters': {
                'n_agents': sim.pars.n_agents,
                'dur':      sim.t.dur,
                'dt':       sim.t.dt,
                },
            'cpu_performance': ratio,
            }

    if verbose:
        if previous:
            print('Previous benchmark:')
            sc.pp(previous)

        print('\nNew benchmark:')
        sc.pp(json)
    else:
        brief = sc.dcp(json['time'])
        brief['cpu_performance'] = json['cpu_performance']
        sc.pp(brief)

    if do_save:
        sc.saveyaml(filename=benchmark_filename, obj=json)

    if verbose:
        print('Done.')

    return json


if __name__ == '__main__':
    do_plot = True
    sc.options(interactive=do_plot)

    T = sc.timer()

    json = test_benchmark(do_save=True) # Run this first so benchmarking is available even if results are different
    new  = test_baseline()
    sim  = make_sim(run=do_plot)

    T.toc()
