from distutils.core import setup
from Cython.Build import cythonize

setup(
    name= 'Fast digestion',
    ext_modules = cythonize("digestfast.pyx"), 
)
