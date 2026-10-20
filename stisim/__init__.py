from .version import __version__, __versiondate__, __license__

from .calibration   import *
from .care_seeking  import *
from .connectors    import *
from .diseases      import *
from .demographics  import *
from .interventions import *
from .logistics     import *
from .networks      import *
from .utils         import *
from .analyzers     import *
from .parameters    import *
from .sim           import *

# Assign the root folder
import sciris as sc
root = sc.thispath(__file__).parent
data = root/'data'

# Double-check key requirements -- should match pyproject.toml
sc.require(['starsim>=3.6.0'], message=f'The following dependencies for STIsim {__version__} were not met: <MISSING>.', die=False)
del sc # Don't keep this in the module

# UPSTREAM: a run sim with declared HIV and ART needs a recursion depth of about 500 to pickle (the expression trees of
# declared modules), which a worker of ss.parallel() inside another deep stack (pytest-xdist, Jupyter) can exceed
import sys as _sys
_sys.setrecursionlimit(max(_sys.getrecursionlimit(), 3000))
del _sys

