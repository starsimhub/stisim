"""
GUD-syphilis connector
"""

from .hiv_sti import connector, effect, where, is_declared

__all__ = ['gud_syph']


@connector(
    rel_sus_syph_gud = 2,   # People with GUD are 2x more likely to acquire syphilis
    rel_trans_syph_gud = 2, # People with GUD are 2x more likely to transmit syphilis
    rel_sus_gud_syph = 2,   # People with primary or secondary syphilis are 2x more likely to acquire GUD
    rel_trans_gud_syph = 2, # People with primary or secondary syphilis are 2x more likely to transmit GUD
)
def gud_syph(gud, syph, p):
    """
    GUD and syphilis: each raises susceptibility to and infectiousness of the other

    Starsim 3's class did nothing (its `update()` was never called); this applies it. Not
    added by `sti.Sim` on its own: pass it as `connectors=`.

    Args:
        gud: the GUD module (`sti.GUD`, or `sti.GUDPlaceholder`, whose `symptomatic` is used)
        syph: the syphilis module
        pars (dict): `rel_sus_syph_gud`, `rel_trans_syph_gud`, `rel_sus_gud_syph`, `rel_trans_gud_syph`

    Examples:
        ```python
        gud, syph = sti.GUD(init_prev=0.05, beta_m2f=0.2), sti.Syphilis()
        sim = sti.Sim(diseases=[gud, syph], connectors=sti.gud_syph(gud, syph))
        ```
    """
    has_gud = where(gud, 'infected' if is_declared(gud) else 'symptomatic')
    active = where(syph, 'symptomatic')
    out = [
        effect('syph_sus_gud', syph, 'rel_sus', p.rel_sus_syph_gud, has_gud),
        effect('syph_trans_gud', syph, 'rel_trans', p.rel_trans_syph_gud, has_gud),
    ]
    if is_declared(gud) or hasattr(gud, 'rel_sus'): # The placeholder has no transmission to modify
        out += [
            effect('gud_sus_syph', gud, 'rel_sus', p.rel_sus_gud_syph, active),
            effect('gud_trans_syph', gud, 'rel_trans', p.rel_trans_gud_syph, active),
        ]
    return out
