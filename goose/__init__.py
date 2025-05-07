"""A Python package for making IDRs and IDR variants"""

# Add imports here
from goose.create import *
from goose.analyze import *
from goose.optimize import *
from goose.backend.optimizer_properties import *
import os 

# Generate _version.py if missing and in the Read the Docs environment
if os.getenv("READTHEDOCS") == "True" and not os.path.isfile('../goose/_version.py'):   
    import versioningit            
    __version__ = versioningit.get_version('../')
else:
    from ._version import __version__

# code that allows access to the data directory
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)

