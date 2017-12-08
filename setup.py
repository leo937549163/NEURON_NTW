from distutils.core import setup
from Cython.Build import cythonize

setup(
	name = 'gperftools_wrapped',
	ext_modules = cythonize("gperftools_wrapped.pyx"),
)
