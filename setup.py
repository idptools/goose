import os
import numpy
from setuptools import Extension, setup, find_packages

# Path to the Cython source and its pre-generated C file.
pyx_file = os.path.join("goose", "backend", "fast_mutations.pyx")
c_file = os.path.join("goose", "backend", "fast_mutations.c")

# Prefer building from the .pyx via Cython when available; fall back to the
# shipped .c file so installs work even if Cython is missing from the build env
# (e.g. --no-build-isolation in conda-forge style builds).
try:
    from Cython.Build import cythonize
    ext_modules = cythonize(
        [
            Extension(
                name="goose.backend.fast_mutations",
                sources=[pyx_file],
                include_dirs=[numpy.get_include()],
            )
        ],
        compiler_directives={"language_level": "3"},
    )
except ImportError:
    if not os.path.exists(c_file):
        raise RuntimeError(
            "Cython is not available and the pre-generated C file "
            f"{c_file} is missing. Install Cython or use a source "
            "distribution that bundles the generated C file."
        )
    ext_modules = [
        Extension(
            name="goose.backend.fast_mutations",
            sources=[c_file],
            include_dirs=[numpy.get_include()],
        )
    ]

setup(
    ext_modules=ext_modules,
    packages=find_packages(),
    include_package_data=True,
)