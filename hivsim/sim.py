import sciris as sc
import starsim as ss
import stisim as sti

__all__ = ['Sim', 'demo']

class Sim(sti.Sim):
    """
    A subclass of sti.Sim that is specifically designed for HIV simulations.

    Currently, this simply parses input parameters among the sim and the HIV module,
    and adds default demographics (pregnancy and deaths), networks (sexual and maternal),
    and interventions (testing, ART, VMMC, and PrEP).

    In future this will support location data and other features.
    """
    def __init__(self, pars=None, sim_pars=None, hiv_pars=None, location=None, **kwargs):

        if location is not None:
            raise NotImplementedError('Location-based sim creation is not implemented yet')

        # Initialize parameters
        pars = sc.mergedicts(pars, kwargs)
        sim_pars = sc.mergedicts(sim_pars)
        hiv_pars = sc.mergedicts(hiv_pars)
        default_sim_keys = ss.SimPars().keys()
        default_hiv_keys = sti.HIVPars().keys()

        # Pull modules out for special processing
        modules = sc.objdict()
        for mod_type in ['diseases', 'networks', 'demographics', 'interventions']:
            modules[mod_type] = sc.mergelists(pars.pop(mod_type, None)) # Remove from kwargs and turn into a list

        # Handle diseases -- first, figure out what parameters belong in HIV
        for key in list(pars.keys()):
            if key in default_hiv_keys:
                if key in default_sim_keys: # If the key is in both, copy
                    val = pars[key]
                else: # Else, pop
                    val = pars.pop(key)
                hiv_pars[key] = val

        hiv = sti.HIV(pars=hiv_pars)
        modules.diseases.insert(0, hiv)

        # Handle demographics
        if not modules.demographics:
            modules.demographics = [ss.Pregnancy(), ss.Deaths()]

        # Handle networks
        if not modules.networks:
            modules.networks = [sti.StructuredSexual(), ss.MaternalNet()]

        # Handle interventions
        if not modules.interventions:
            modules.interventions = [sti.HIVTest(), sti.ART(), sti.VMMC(), sti.Prep()]

        super().__init__(
            pars          = sim_pars,
            demographics  = modules.demographics,
            networks      = modules.networks,
            diseases      = modules.diseases,
            interventions = modules.interventions,
            **pars
        )
        return


# Available examples for hivsim.demo()
EXAMPLES = {
    'simple':   'Minimal HIV sim with hivsim defaults',
    'zimbabwe': 'Zimbabwe HIV model with calibrated parameters and UNAIDS data',
}


def demo(example=None, run=True, plot=True, **kwargs):
    """
    Create a demo HIVsim simulation.

    Args:
        example (str): Example name ('simple', 'zimbabwe'). Default: 'simple'.
        run (bool): Whether to run the sim.
        plot (bool): Whether to plot results (only if run=True).
        **kwargs: Passed to the example's make_sim().

    Returns:
        Sim: Configured (and optionally run) simulation.

    Examples::

        import hivsim
        hivsim.demo()                                     # Run simple default demo
        hivsim.demo('zimbabwe')                            # Run Zimbabwe HIV model
        sim = hivsim.demo('zimbabwe', run=False, n_agents=500)  # Just create it
    """
    if example is None:
        example = 'simple'
    if example not in EXAMPLES:
        available = ', '.join(EXAMPLES.keys())
        raise ValueError(f"Example '{example}' not found. Available: {available}")

    mod = sc.importbyname(f'hivsim_examples.{example}.sim')
    sim = mod.make_sim(**kwargs)

    if run:
        sim.run()
        if plot:
            try:
                sim.plot('hiv', annualize=True)
            except Exception:
                sim.plot('hiv')

    return sim