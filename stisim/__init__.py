from .version import __version__, __versiondate__, __license__

from .calibration   import *
from .care_seeking  import *
from .connectors    import *
from .diseases      import *
from .demographics  import *
from .interventions import *
from .networks      import *
from .utils         import *
from .analyzers     import *
from .parameters    import *
from .sim           import *
from . import hivsim # Do not import subfunctions

# Assign the root folder
import sciris as sc
root = sc.thispath(__file__).parent
data = root/'data'

# Double-check key requirements -- should match pyproject.toml
sc.require(['starsim>=3.0.0'], message=f'The following dependencies for STIsim {__version__} were not met: <MISSING>.', die=False)
del sc # Don't keep this in the module

