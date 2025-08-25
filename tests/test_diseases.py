"""
Test epi dynamics
"""
import sciris as sc
import starsim as ss
import stisim as sti

n_agents = 2000
debug = False


def test_syph_epi():
    sc.heading('Test epi dynamics of syphilis')

    base_pars = dict(n_agents=n_agents)

    # Define the parameters to vary
    par_effects = dict(
        rel_trans_primary=[0.5, 1.5],  # Relative transmissibility during primary infection
        init_prev=[0.1, 0.9],
        beta_m2f=[0.2, 0.5]  # Beta for male to female transmission; opposite direction uses half this value
    )

    # Loop over each of the above parameters and make sure they affect the epi dynamics in the expected ways
    for par, par_val in par_effects.items():
        lo = par_val[0]
        hi = par_val[1]

        # Make baseline pars
        pars0 = sc.dcp(base_pars)
        pars1 = sc.dcp(base_pars)

        simpardict_lo = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: lo})
        simpardict_hi = sc.mergedicts(dict(beta_m2f=0.3, init_prev=0.05), {par: hi})

        pars0['diseases'] = sti.Syphilis(**simpardict_lo)
        pars1['diseases'] = sti.Syphilis(**simpardict_hi)

        # Run the simulations and pull out the results
        s0 = sti.Sim(pars0, label=f'{par} {par_val[0]}')
        s1 = sti.Sim(pars1, label=f'{par} {par_val[1]}')
        ss.parallel(s0, s1, debug=debug)

        # Check results
        ind = 1 if par == 'init_prev' else -1
        v0 = s0.results.syph.cum_infections[ind]
        v1 = s1.results.syph.cum_infections[ind]

        print(f'Checking with varying {par:10s} ... ', end='')
        assert v0 <= v1, f'Expected infections to be lower with {par}={lo} than with {par}={hi}, but {v0} > {v1})'
        print(f'✓ ({v0} <= {v1})')

    return s0, s1


def test_hiv_epi():
    sc.heading('Test epi dynamics of hiv')

    base_pars = dict(n_agents=n_agents, diseases='hiv', beta_m2f=0.05, init_prev=0.05)

    # Define the parameters to vary
    par_effects = dict(
        # dur_acute=[1, 24],
        # init_prev=[0.01, 0.1],
        beta_m2f=[0.01, 0.2],
    )

    # Loop over each of the above parameters and make sure they affect the epi dynamics in the expected ways
    for par, par_val in par_effects.items():
        lo = par_val[0]
        hi = par_val[1]

        # Make baseline pars
        pars0 = sc.dcp(base_pars)
        pars1 = sc.dcp(base_pars)

        pars0[par] = lo
        pars1[par] = hi

        # Run the simulations and pull out the results
        s0 = sti.Sim(pars0, label=f'{par} {par_val[0]}')
        s1 = sti.Sim(pars1, label=f'{par} {par_val[1]}')
        ss.parallel(s0, s1, debug=debug)

        # Check results
        ind = 1 if par == 'init_prev' else -1
        v0 = s0.results.hiv.cum_infections[ind]
        v1 = s1.results.hiv.cum_infections[ind]

        print(f'Checking with varying {par:10s} ... ', end='')
        assert v0 <= v1, f'Expected infections to be lower with {par}={lo} than with {par}={hi}, but {v0} > {v1})'
        print(f'✓ ({v0} <= {v1})')

    return s0, s1


if __name__ == '__main__':
    sc.options(interactive=False)
    s1, s2 = test_syph_epi()
    s3, s4 = test_hiv_epi()
