
Welcome to Mrsimulator's documentation!
=======================================

.. only:: html

    .. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master
        :target: https://travis-ci.org/DeepanshS/mrsimulator
        :alt: Build Status

    .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
        :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
        :alt: Documentation Status

    .. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
        :target: https://opensource.org/licenses/BSD-3-Clause
        :alt: License

    .. image:: https://img.shields.io/github/issues-raw/deepanshs/mrsimulator
        :target: https://github.com/DeepanshS/mrsimulator/issues
        :alt: GitHub issues

    .. image:: https://badge.fury.io/py/mrsimulator.svg
        :target: https://badge.fury.io/py/mrsimulator
        :alt: PyPI version

    .. image:: https://img.shields.io/pypi/dm/mrsimulator.svg
        :target: https://img.shields.io/pypi/dm/mrsimulator
        :alt: PyPI - Downloads

    .. image:: https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18
        :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/
        :alt: Total alerts

    .. image:: https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18
        :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python
        :alt: Language grade: Python

    .. image:: https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/DeepanshS/mrsimulator

----

``Mrsimulator`` is a python package with methods and tools for simulating
fast solid-state nuclear magnetic resonance (NMR) line-shapes. The
core methods is written in C, wrapped and made available
in python for python users.

.. warning::
    The package is currently under development. We advice using with
    caution. Bug report are greatly appreciated.

----

**Features**

At present, ``mrsimulator`` offers fast-simulation of one-dimensional NMR
line-shape of uncoupled spins for the following:

- Spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}`,
  (See the list of supported isotopes),
- Arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin frequency,
- Variable angle spinning (VAS) at arbitrary angle and spin frequency,
- Static line-shape.

----

    **View the example gallery**

    .. only:: html

        .. image:: https://img.shields.io/badge/View-Example%20Gallery-Purple?size=large
            :target: auto_examples/index.html

----

**Contribution**

``Mrsimulator`` is a open source NMR simulation package. We are a small team
working on developing the package for the NMR community. Any contribution and
suggestion is greatly appreciated.


-------

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   installation
   requirements
   understanding-isotopomers
   getting_started
   using_mrsimulator_objects
   load_isotopomers
   configuring_simulator
   benchmark
   auto_examples/index
   theory/components
   api_py/py_api
   api_c/c_api


..    objects
..    spectrum_object
..    theory/wigner_rotations
..    examples


.. only:: html

    Indices and tables
    ------------------

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
