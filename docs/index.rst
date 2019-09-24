
.. only:: html

    .. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master
        :target: https://travis-ci.org/DeepanshS/mrsimulator
    .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
        :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
        :alt: Documentation Status

Welcome to mrsimulator's documentation!
=======================================

The package ``mrsimulator`` is a collection of methods and tools for simulating
nuclear magnetic resonance (NMR) line-shapes. Because the bulk of the code is
written in C and eventually wrapped in python, ``mrsimulator`` is a fast NMR
line-shape simulation library, harnessing both the speed of C and the
advantages of Python programming language.

The package is currently in development. We currently support simulation of
NMR line-shapes from uncoupled spin :math:`I=\frac{1}{2}` nucleus under static,
magic angle spinning (MAS), and variable angle spinning (VSA) conditions.

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
