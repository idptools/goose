import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup, find_packages
import os

# set cython file
cython_file = os.path.join("goose", "backend", "fast_mutations.pyx")

# Define the extension
extensions = [
    Extension(
        name="goose.backend.fast_mutations",
        sources=[cython_file],
        include_dirs=[numpy.get_include()],
    )
]

# Setup function
setup(
    ext_modules = cythonize(extensions, compiler_directives={'language_level' : "3"}),
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True
    )