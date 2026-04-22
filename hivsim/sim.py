import sciris as sc
import starsim as ss
import stisim as sti

__all__ = ['Sim', 'demo']


# Defaults supplied by hivsim when the user doesn't override them.
# Exposed at module level so callers can see exactly what they're getting.
def _default_demographics():
    return [ss.Pregnancy(), ss.Deaths()]


def _default_interventions():
    return [sti.HIVTest(), sti.ART(), sti.VMMC(), sti.Prep()]


def _default_networks(demographics):
    """Default networks: sexual + maternal, + breastfeeding if Pregnancy is present."""
    nets = [sti.StructuredSexual(), ss.MaternalNet()]
    if any(isinstance(m, ss.Pregnancy) for m in sc.tolist(demographics)):
        nets.append(ss.BreastfeedingNet())
    return nets


class Sim(sti.Sim):
    """
    A thin subclass of sti.Sim that supplies HIV-appropriate defaults.

    Equivalent to::

        sti.Sim(
            diseases='hiv',
            demographics=[ss.Pregnancy(), ss.Deaths()],
            networks=[sti.StructuredSexual(), ss.MaternalNet(), ss.BreastfeedingNet()],
            interventions=[sti.HIVTest(), sti.ART(), sti.VMMC(), sti.Prep()],
        )

    Parameter routing (hiv_pars, network pars, etc.) is entirely delegated to
    sti.Sim's ``separate_pars`` machinery. Users can customize:
        - A whole module slot by passing the module directly (``diseases=MyHIV()``)
        - Or by passing pars for the default module (``hiv=dict(rel_trans_acute=26)``
          or flat ``rel_trans_acute=26``)
        - But not both for the same slot (sti.Sim raises XOR violation).

    Collisions between ``pars`` and ``kwargs`` raise an error: specify a value
    exactly once.
    """
    def __init__(self, pars=None, sim_pars=None, hiv_pars=None, location=None, **kwargs):

        if location is not None:
            raise NotImplementedError('Location-based sim creation is not implemented yet')

        # Fail loudly on overlapping keys between `pars` and `**kwargs` — rather
        # than silently taking the later value, force the caller to pick one.
        if pars is not None:
            overlap = set(pars) & set(kwargs)
            if overlap:
                raise ValueError(
                    f'Keys appear in both `pars` and kwargs: {sorted(overlap)}. '
                    f'Specify each parameter in exactly one place.'
                )

        pars = sc.mergedicts(pars, kwargs)

        # Default diseases → 'hiv' string, so sti.Sim's string-path constructs
        # HIV and applies routed sti_pars (including `hiv=dict(...)` or flat keys).
        pars.setdefault('diseases', 'hiv')

        # Default demographics → pregnancy + deaths.
        demographics = pars.pop('demographics', None) or _default_demographics()

        # Default networks — since there's no string shortcut for our combo,
        # consume network pars from kwargs into the instance here. (If the user
        # passed their own `networks=[...]`, any nw_pars they also passed will
        # still be routed by sti.Sim and caught by its XOR validation.)
        if 'networks' not in pars:
            nw_kwargs = dict(pars.pop('nw_pars', None) or {})
            for k in list(pars):
                if k in sti.NetworkPars().keys():
                    nw_kwargs[k] = pars.pop(k)
            if nw_kwargs:
                # Apply user's nw_pars to the default StructuredSexual.
                pars['networks'] = [sti.StructuredSexual(**nw_kwargs), ss.MaternalNet()]
                if any(isinstance(m, ss.Pregnancy) for m in sc.tolist(demographics)):
                    pars['networks'].append(ss.BreastfeedingNet())
            else:
                pars['networks'] = _default_networks(demographics)

        # Default interventions.
        pars.setdefault('interventions', _default_interventions())

        # Backward-compat: hiv_pars kwarg → merge into the 'hiv' keyed dict
        # so sti.Sim's separate_pars picks it up.
        if hiv_pars:
            pars['hiv'] = sc.mergedicts(pars.get('hiv'), hiv_pars)

        super().__init__(pars=sim_pars, demographics=demographics, **pars)
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
