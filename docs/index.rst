
.. only:: html

    .. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master
        :target: https://travis-ci.org/DeepanshS/mrsimulator
    .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
        :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
        :alt: Documentation Status

Welcome to mrsimulator's documentation!
=======================================

The package, ``mrsimulator``, is a library of methods and tools for fast
simulation of solid-state nuclear magnetic resonance (NMR) line-shapes.
The library is written in C and is wrapped and made available in python for
python users. The package is currently under development

Features
********

At present, ``mrsimulator`` offers fast-simulation of one-dimensional NMR
line-shape of uncoupled sites for the following:

- Spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}`,
  (See the list of supported isotopes),
- Arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin frequency,
- Variable angle spinning (VAS) at arbitrary angle and spin frequency,
- Static line-shape.

-------

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
