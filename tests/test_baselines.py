"""
Test that the current version of HPVsim exactly matches
the baseline results.
"""

import numpy as np
import sciris as sc
import stisim as sti
import starsim as ss

do_save = 1
baseline_filename  = sc.thisdir(__file__, 'baseline.json')
benchmark_filename = sc.thisdir(__file__, 'benchmark.json')


def make_sim(**kwargs):
    '''
    Define a default simulation for testing the baseline, including
    interventions to increase coverage. If run directly (not via pytest), also
    plot the sim by default.
    '''
    bv = sti.BV()
    ct = sti.Chlamydia()
    gon = sti.Gonorrhea()
    gud = sti.GUD()
    hiv = sti.HIV()
    st = sti.Syphilis()
    trich = sti.Trichomoniasis()

    diseases = [bv, ct, gon, gud, hiv, st, trich]
    networks = [sti.FastStructuredSexual()]
    connectors = [sti.gud_syph(), sti.hiv_syph(hiv, st), sti.hiv_bv(hiv, bv), sti.hiv_ct(hiv, ct), sti.hiv_tv(hiv, trich)]
    interventions = [sti.ART(), sti.SyphTx()]
    demographics = [sti.Pregnancy(), ss.Deaths()]



    # Define the parameters
    pars = dict(
        n_agents      = 10e3,       # Population size
        start         = 2000,       # Starting year
        stop          = 2010,         # Number of years to simulate
        verbose       = 0,          # Don't print details of the run
        rand_seed     = 2,          # Set a non-default seed
    )

    # Create the sim
    sim = ss.Sim(pars=pars, diseases=diseases, networks=networks, connectors=connectors, interventions=interventions, demographics=demographics)

    return sim


def save_baseline():
    '''
    Refresh the baseline results. This function is not called during standard testing,
    but instead is called by the update_baseline script.
    '''

    print('Updating baseline values...')

    # Export results
    s = make_sim()
    s.run()
    s.to_json(filename=baseline_filename)

    print('Done.')

    return


def test_baseline():
    ''' Compare the current default sim against the saved baseline '''

    # Load existing baseline
    baseline = sc.loadjson(baseline_filename)
    old = baseline['summary']

    # Calculate new baseline
    new = make_sim()
    new.run()

    # Compute the comparison
    ss.diff_sims(old, new, full=True, die=True)

    return new


def test_benchmark(do_save=do_save, repeats=1, verbose=True):
    ''' Compare benchmark performance '''

    if verbose: print('Running benchmark...')
    previous = sc.loadjson(benchmark_filename)

    t_inits = []
    t_runs  = []

    def normalize_performance():
        ''' Normalize performance across CPUs '''
        t_bls = []
        bl_repeats = 3
        n_outer = 10
        n_inner = 1e6
        for r in range(bl_repeats):
            t0 = sc.tic()
            for i in range(n_outer):
                a = np.random.random(int(n_inner))
                b = np.random.random(int(n_inner))
                a*b
            t_bl = sc.toc(t0, output=True)
            t_bls.append(t_bl)
        t_bl = min(t_bls)
        reference = 0.07 # Benchmarked on an Intel i7-12700H CPU @ 2.90GHz
        ratio = reference/t_bl
        return ratio


    # Test CPU performance before the run
    r1 = normalize_performance()

    # Do the actual benchmarking
    for r in range(repeats):

        # Create the sim
        sim = make_sim(verbose=0)

        # Time initialization
        t0 = sc.tic()
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
    r2 = normalize_performance()
    ratio = (r1+r2)/2
    t_init = min(t_inits)*ratio
    t_run  = min(t_runs)*ratio

    # Construct json
    n_decimals = 3
    json = {'time': {
                'initialize': round(t_init, n_decimals),
                'run':        round(t_run,  n_decimals),
                },
            'parameters': {
                'n_agents': sim.pars['n_agents'],
                'n_years':   sim.pars['stop']-sim.pars['start'],
                'n_diseases':  len(sim.diseases),
                'n_interventions': len(sim['interventions']),
                'n_analyzers': len(sim['analyzers']),
                },
            'cpu_performance': ratio,
            }

    if verbose:
        print('Previous benchmark:')
        sc.pp(previous)

        print('\nNew benchmark:')
        sc.pp(json)
    else:
        brief = sc.dcp(json['time'])
        brief['cpu_performance'] = json['cpu_performance']
        sc.pp(brief)

    if do_save:
        sc.savejson(filename=benchmark_filename, obj=json, indent=2)

    if verbose:
        print('Done.')

    return json



if __name__ == '__main__':

    T = sc.tic()

    json = test_benchmark(do_save=do_save, repeats=5, verbose=False) # Run this first so benchmarking is available even if results are different
    new  = test_baseline()
    make_sim()

    print('\n'*2)
    sc.toc(T)
    print('Done.')