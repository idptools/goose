"""A Python package for making IDRs and IDR variants"""

# Add imports here
from .create import *
from .analyze import *
from .optimize import *
from .backend.optimizer_properties import *
import os 

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

# code that allows access to the data directory
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)

