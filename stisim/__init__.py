# Last verified on: Starsim 4.0.0rc4 (34266c4c4)

from .version import __version__, __versiondate__, __license__

from .utils         import *
from .parameters    import *
from .demographics  import *
from .networks      import *
from .diseases      import *
from .connectors    import *
from .care_seeking  import *
from .logistics     import *
from .interventions import *
from .analyzers     import *
from .sim           import *

# Assign the root and data folders
import sciris as sc
root = sc.thispath(__file__).parent
data = root/'data'
del sc # Don't keep this in the module
