"""STIsim networks package.

See ``README.md`` for the layout. Public API exports are flat: e.g.
``stisim.MFNetwork`` works, as does ``stisim.networks.MFNetwork``.

For direct access to matcher functions used by ``MFNetwork.match_pairs``,
use ``stisim.networks.matchers`` (see Commit B).
"""
from .base import NoPartnersFound, BasePars, NetworkPars, BaseNetwork
from .mf import MFPars, MFNetwork
from .fsw import SWPars, SWNetwork
from .msm import AgeMatchedMSM, AgeApproxMSM
from .layered_networks import StructuredSexual, PriorPartners
from .scalefree import MSMScaleFreeNetwork
from . import matchers
from .matchers import MATCHERS

__all__ = [
    'NoPartnersFound', 'BasePars', 'NetworkPars', 'BaseNetwork',
    'MFPars', 'MFNetwork',
    'SWPars', 'SWNetwork',
    'AgeMatchedMSM', 'AgeApproxMSM',
    'StructuredSexual', 'PriorPartners',
    'MSMScaleFreeNetwork',
    'MATCHERS',
]
