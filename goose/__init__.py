"""A Python package for making IDRs and IDR variants"""

# Add imports here
from .create import *
from .analyze import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
