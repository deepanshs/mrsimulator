
.. only:: html

    .. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master
        :target: https://travis-ci.org/DeepanshS/mrsimulator
    .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
        :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
        :alt: Documentation Status

Welcome to mrsimulator's documentation!
=======================================

``mrsimulator`` is a library package with methods and tools for fast
simulation of solid-state nuclear magnetic resonance (NMR) line-shapes.
The library contains routines written in C which are wrapped and made
available in python.

The package is currently under development. At present, `mrsimulator` features
simulation of one-dimensional NMR line-shape of uncoupled spin
:math:`I=\frac{1}{2}` isotopes for the following scenarios --

- At arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin rate,
- Variable angle spinning (VAS) at arbitrary angle and spin rates,
- Static line-shape.

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   installation
   requirements
   getting_started
   load_isotopomers
   objects
   spectrum_object
   examples
   theory/components
   api_py/py_api
   api_c/c_api


..    theory/wigner_rotations


.. only:: html

    Indices and tables
    ------------------

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
