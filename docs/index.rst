.. MRSimulator documentation master file, created by
   sphinx-quickstart on Tue Apr  2 21:31:08 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master
    :target: https://travis-ci.org/DeepanshS/mrsimulator
.. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
    :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status

Welcome to mrsimulator's documentation!
=======================================

The package ``mrsimulator`` is a collection of methods and tools for simulating
nuclear magnetic resonance (NMR) lineshapes. Because the bulk of the code is
written in C and eventually wrapped in python, ``mrsimulator`` is a fast NMR
lineshape simulation library, harnessing both the speed of C and the advantages
of Python programming language.

The package is currently in development. We currently support simulation of
NMR lineshapes from uncoupled spin :math:`I=\frac{1}{2}` nucleus under static,
magic angle spinning (MAS), and variable angle spinning (VSA) conditions.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   mr_objects
   installation
   getting_started
   load_isotopomers
   examples
   what's_going_on
   api_py/py_api
   api_c/c_api
   c_api/c_api

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
