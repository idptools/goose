import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

# Define the extension
extensions = [
    Extension(
        "goose.backend.fast_mutations",
        sources=["goose/backend/fast_mutations.pyx"],
        include_dirs=[numpy.get_include()],
    )
]

# Setup function
setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": 3,
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
    ),
)
