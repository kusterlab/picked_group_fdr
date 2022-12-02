# See if Cython is installed
try:
    from Cython.Build import cythonize
# Do nothing if Cython is not available
except ImportError:
    # Got to provide this function. Otherwise, poetry will fail
    def build(setup_kwargs):
        pass
# Cython is installed. Compile
else:
    from distutils.core import setup
    
    import numpy
    
    # This function will be executed in setup.py:
    def build(setup_kwargs):
        setup(
            name= 'PickedGroupFDR',
            ext_modules = cythonize(["picked_group_fdr/digestfast.pyx", "picked_group_fdr/quant/lfq_helpers.pyx"]), 
            include_dirs=[numpy.get_include()],
        )
