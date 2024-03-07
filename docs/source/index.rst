.. Picked group FDR documentation master file, created by
   sphinx-quickstart on Wed Feb 21 15:34:00 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Picked Group FDR's documentation!
============================================


|PyPI version| |Supported Python versions| |PyPI downloads|

Scalable, accurate and sensitive protein group FDRs for large-scale mass
spectrometry experiments.

Different search engine outputs are supported:

* `MaxQuant <https://www.maxquant.org/>`__ (LFQ, TMT, SILAC)
* `Percolator <https://github.com/percolator/percolator>`__ (no quantification)
* `FragPipe <https://fragpipe.nesvilab.org/>`__ (LFQ)
* `Sage <https://github.com/lazear/sage>`__ (LFQ)

Examples on how to run Picked Protein Group FDR with these different
search engine outputs can be found under :doc:`Basic usage <usage>` and in the 
``data/`` folder.


.. |PyPI version| image:: https://img.shields.io/pypi/v/picked_group_fdr.svg?logo=pypi&logoColor=FFE873
   :target: https://pypi.org/project/picked_group_fdr/
.. |Supported Python versions| image:: https://img.shields.io/pypi/pyversions/picked_group_fdr.svg?logo=python&logoColor=FFE873
   :target: https://pypi.org/project/picked_group_fdr/
.. |PyPI downloads| image:: https://img.shields.io/pypi/dm/picked_group_fdr.svg
   :target: https://pypistats.org/packages/picked_group_fdr


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage
   manuscript
   cite