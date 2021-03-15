import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name='CYGNO reconstruction',
  ext_modules=[Extension('_ddbscan_inner_cy', ['cluster/ddbscan_inner_cython.pyx'],)],
  cmdclass={'build_ext': build_ext},
)
