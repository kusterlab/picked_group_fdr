from distutils.core import setup
from Cython.Build import cythonize

import numpy 

setup(
    name= 'PickedGroupFDR',
    ext_modules = cythonize(["picked_group_fdr/digestfast.pyx", "picked_group_fdr/quant/lfq_helpers.pyx"]), 
    include_dirs=[numpy.get_include()],
)
