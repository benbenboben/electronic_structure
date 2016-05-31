

import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np                           # <---- New line
import scipy as la                           # <---- New line
import copy

#os.environ["CC"] = "clang"
#os.environ["CXX"] = "clang++"
os.environ["CC"] = "gcc-5"
os.environ["CXX"] = "g++-5"
#os.environ["CC"] = "icc"
#os.environ["CXX"] = "icpc"

# for GNU
os.environ["ARCHFLAGS"] = "-arch x86_64"

ext_modules = [Extension("ao2mo", ["ao2mo.pyx"], extra_compile_args=['-fopenmp', '-O3'], extra_link_args=['-fopenmp'])]

setup(
  name = 'ao2mo.pyx',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include(), la.get_include()],         # <---- New line
  ext_modules = ext_modules
  )
