"""
Coinfection connectors, as lists of modifiers

Each connector is a function of two diseases that returns a list of `ss.modify()`s: one for
each effect of one disease on the other's susceptibility or infectiousness. Pass the list (or
several, joined) as `connectors=`; `sti.Sim` adds the connector for each pair of its diseases
that has one (see `sti.connector_register()`).
"""

import functools
import numpy as np
import sciris as sc
import starsim as ss

__all__ = ['connector', 'legacy_modify', 'hiv_syph', 'hiv_tv', 'hiv_ng', 'hiv_ct', 'hiv_bv', 'hiv_simplebv']


def is_declared(mod):
    return getattr(mod, 'graph', None) is not None


class where:
    """
    Who has a state of a disease (e.g. `where(syph, 'symptomatic')`): a lazy condition if both diseases are declared,
    else a function of the sim that reads the module's array of that name as agents (for a Starsim 3 module)
    """
    def __init__(self, mod, attr):
        self.mod = mod
        self.attr = attr
        return

    def cond(self, target):
        if is_declared(self.mod) and is_declared(target):
            return getattr(self.mod, self.attr)
        name, attr = self.mod.name, self.attr
        def has(sim):
            return np.asarray(getattr(sim.diseases[name], attr).raw, dtype=bool)
        has.__name__ = f'{name}_{attr}'
        return has


def cd4_below(hiv, threshold):
    """ People with HIV whose CD4 count is below a threshold, as a function of the sim. TODO: an attribute condition (hiv.cd4 < 500) once HIV is declared """
    name = hiv.name
    def below(sim):
        cd4 = sim.diseases[name].cd4.raw
        with np.errstate(invalid='ignore'):
            return np.asarray(cd4 < threshold) & ~np.isnan(cd4)
    below.__name__ = f'cd4_below_{threshold}'
    return below


class legacy_modify(ss.Connector):
    """
    Multiply `rel_sus` or `rel_trans` of a Starsim 3 module where a condition holds, as `ss.modify()` does for a
    declared one; the module resets them at each step (in `step_state()`). TODO(v4-rs): ss.modify() of a Starsim 3 module

    Args:
        target (str): the name of the module
        kind (str): `'rel_sus'` or `'rel_trans'`
        multiply (float): the multiplier
        where (func): a function of the sim returning who it applies to (a boolean array over all agents, or UIDs)
        name (str): the name of the connector
    """
    def __init__(self, target, kind, multiply, where=None, name=None, **kwargs):
        super().__init__(name=name, **kwargs)
        self.target_name = getattr(target, 'name', target)
        self.kind = kind
        self.multiply = multiply
        self.where = where
        return

    def step(self):
        arr = getattr(self.sim.diseases[self.target_name], self.kind)
        if self.where is None:
            arr[:] *= self.multiply
            return
        sel = self.where(self.sim)
        uids = ss.uids(sel) if isinstance(sel, ss.uids) else ss.uids(np.flatnonzero(np.asarray(sel, dtype=bool)))
        uids = uids[self.sim.people.alive.raw[uids]]
        arr[uids] = arr[uids]*self.multiply
        return


def effect(name, target, kind, multiply, cond):
    """ One effect of a coinfection: `ss.modify()` of a declared disease's rel_sus or rel_trans, or `sti.legacy_modify` of a Starsim 3 one """
    if isinstance(cond, where):
        cond = cond.cond(target)
    if is_declared(target):
        return ss.modify(getattr(target, kind), multiply=multiply, where=cond, name=name)
    return legacy_modify(target, kind, multiply=multiply, where=cond, name=name)


def connector(**defaults):
    """
    Make a coinfection connector from a function `fn(d1, d2, pars)` that returns a list of modifiers

    The connector is called as `conn(d1, d2, pars=None, **kwargs)`, with parameters by name (the defaults are
    `conn.pars`), and returns the list. Calling it with `None` for the diseases returns `conn.pars`, as Starsim 3's
    connector classes did for `sti.Sim`'s parameter routing.
    """
    def decorator(fn):
        @functools.wraps(fn)
        def make(d1=None, d2=None, pars=None, **kwargs):
            given = sc.mergedicts(pars, kwargs)
            unknown = set(given) - set(defaults)
            if unknown:
                errormsg = f'Unknown parameters for {fn.__name__}: {sc.strjoin(sorted(unknown))}; the parameters are {sc.strjoin(defaults.keys())}'
                raise ValueError(errormsg)
            p = sc.objdict(sc.mergedicts(defaults, given))
            if d1 is None and d2 is None:
                return p
            return [mod for mod in fn(d1, d2, p) if mod is not None]
        make.pars = sc.objdict(defaults)
        make.is_connector = True
        return make
    return decorator


@connector(
    rel_sus_hiv_syph = 2.67,  # Relative susceptibility to HIV with primary or secondary syphilis
    rel_trans_hiv_syph = 1.2, # Relative HIV infectiousness with primary or secondary syphilis
    rel_sus_syph_hiv = 1,     # Relative susceptibility to syphilis with HIV (CD4 < 500)
    rel_sus_syph_aids = 1,    # Relative susceptibility to syphilis with AIDS (CD4 < 200)
    rel_trans_syph_hiv = 1,   # Relative syphilis infectiousness with HIV (CD4 < 500)
    rel_trans_syph_aids = 1,  # Relative syphilis infectiousness with AIDS (CD4 < 200)
)
def hiv_syph(hiv, syph, p):
    """
    HIV and syphilis: primary or secondary syphilis raises susceptibility to HIV and HIV
    infectiousness; low CD4 can raise susceptibility to and infectiousness of syphilis

    Args:
        hiv: the HIV module
        syph: the syphilis module (or `sti.SyphilisPlaceholder`, which only has the effects on HIV)
        pars (dict): `rel_sus_hiv_syph`, `rel_trans_hiv_syph`, `rel_sus_syph_hiv`, `rel_sus_syph_aids`, `rel_trans_syph_hiv`, `rel_trans_syph_aids`

    Examples:
        ```python
        hiv, syph = sti.HIV(), sti.Syphilis()
        sim = sti.Sim(diseases=[hiv, syph], connectors=sti.hiv_syph(hiv, syph, rel_sus_hiv_syph=3))
        ```
    """
    active = where(syph, 'symptomatic')
    out = [
        effect('hiv_sus_syph', hiv, 'rel_sus', p.rel_sus_hiv_syph, active),
        effect('hiv_trans_syph', hiv, 'rel_trans', p.rel_trans_hiv_syph, active),
    ]
    if hasattr(syph, 'active') and not is_declared(syph): # The placeholder
        return out
    low, aids = cd4_below(hiv, 500), cd4_below(hiv, 200)
    out += [
        effect('syph_sus_hiv', syph, 'rel_sus', p.rel_sus_syph_hiv, low),
        effect('syph_sus_aids', syph, 'rel_sus', p.rel_sus_syph_aids, aids),
        effect('syph_trans_hiv', syph, 'rel_trans', p.rel_trans_syph_hiv, low),
        effect('syph_trans_aids', syph, 'rel_trans', p.rel_trans_syph_aids, aids),
    ]
    return out


@connector(
    rel_sus_hiv_tv = 1.5, # TV: 1.5x risk of HIV acquisition (https://www.who.int/news-room/fact-sheets/detail/trichomoniasis, https://pubmed.ncbi.nlm.nih.gov/30341233/)
)
def hiv_tv(hiv, tv, p):
    """
    HIV and trichomoniasis: an established TV infection (not the latent period) raises susceptibility to HIV

    Args:
        hiv: the HIV module
        tv: the trichomoniasis module
        pars (dict): `rel_sus_hiv_tv`
    """
    return [effect('hiv_sus_tv', hiv, 'rel_sus', p.rel_sus_hiv_tv, where(tv, 'infectious'))]


@connector(
    rel_sus_hiv_ng = 1.2,   # Having NG increases risk of HIV acquisition
    rel_trans_hiv_ng = 1.2, # Having NG increases HIV transmissibility
)
def hiv_ng(hiv, ng, p):
    """
    HIV and gonorrhea: an established NG infection raises susceptibility to HIV and HIV infectiousness

    Args:
        hiv: the HIV module
        ng: the gonorrhea module
        pars (dict): `rel_sus_hiv_ng`, `rel_trans_hiv_ng`
    """
    ng_inf = where(ng, 'infectious')
    return [
        effect('hiv_sus_ng', hiv, 'rel_sus', p.rel_sus_hiv_ng, ng_inf),
        effect('hiv_trans_ng', hiv, 'rel_trans', p.rel_trans_hiv_ng, ng_inf),
    ]


@connector(
    rel_sus_hiv_ct = 1, # Placeholder: no effect
)
def hiv_ct(hiv, ct, p):
    """
    HIV and chlamydia: an established CT infection changes susceptibility to HIV (no effect by default)

    Args:
        hiv: the HIV module
        ct: the chlamydia module
        pars (dict): `rel_sus_hiv_ct`
    """
    return [effect('hiv_sus_ct', hiv, 'rel_sus', p.rel_sus_hiv_ct, where(ct, 'infectious'))]


@connector(
    rel_sus_hiv_bv = 2,   # Placeholder: BV leads to 2x risk of HIV acquisition
    rel_trans_hiv_bv = 2, # Placeholder: BV leads to 2x risk of HIV transmission
)
def hiv_bv(hiv, bv, p):
    """
    HIV and BV: BV (CST 4) raises susceptibility to HIV and HIV infectiousness

    Args:
        hiv: the HIV module
        bv: the BV module (`sti.BV`; or `sti.SimpleBV`, whose `infected` is used)
        pars (dict): `rel_sus_hiv_bv`, `rel_trans_hiv_bv`
    """
    cst4 = where(bv, 'cst4' if is_declared(bv) else 'infected')
    return [
        effect('hiv_sus_bv', hiv, 'rel_sus', p.rel_sus_hiv_bv, cst4),
        effect('hiv_trans_bv', hiv, 'rel_trans', p.rel_trans_hiv_bv, cst4),
    ]


hiv_simplebv = hiv_bv # Starsim 3 name
