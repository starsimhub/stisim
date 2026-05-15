"""
Minimal HIV simulation example using the heterosexual MFNetwork only.

Use this as a starting point for HIV models that don't need to model
sex-work transmission as a separate compartment. PrEP (which defaults to
FSW eligibility) is not included by default here — add it explicitly with
your own eligibility rule if needed.
"""
import sciris as sc
import starsim as ss
import stisim as sti
import hivsim


def make_sim(**kwargs):
    """Create a simple HIV simulation using ``sti.MFNetwork``.

    Flat ``**kwargs`` are auto-routed via ``sti.route_pars`` — pass network
    pars (``debut_f=18``, ``prop_f0=0.7``), sim pars (``n_agents=5000``,
    ``dur=10``), or sti pars (``init_prev=0.05``) directly.
    """
    kwargs.setdefault('n_agents', 2000)
    kwargs.setdefault('dur', 20)
    nw_kwargs = sti.route_pars(kwargs).nw
    for k in nw_kwargs: kwargs.pop(k)  # consumed here, don't re-route in sti.Sim

    if 'networks' not in kwargs:
        nets = [sti.MFNetwork(**nw_kwargs), ss.MaternalNet()]
        # BreastfeedingNet only makes sense if Pregnancy is in demographics.
        # User-supplied demographics: respect them. Otherwise hivsim.Sim defaults
        # to Pregnancy + Deaths, so include BreastfeedingNet.
        demog = kwargs.get('demographics')
        if demog is None or any(isinstance(m, ss.Pregnancy) for m in sc.tolist(demog)):
            nets.append(ss.BreastfeedingNet())
        kwargs['networks'] = nets
    kwargs.setdefault('interventions', [sti.HIVTest(), sti.ART(), sti.VMMC()])
    return hivsim.Sim(**kwargs)
