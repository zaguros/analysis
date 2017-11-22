from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

import os
print os.getcwd()
setup(
    ext_modules=cythonize("pop_montecarlo_c.pyx"),
    include_dirs=[numpy.get_include()]
)    