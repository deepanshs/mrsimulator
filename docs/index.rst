
================================================
Welcome to the Mrsimulator project documentation
================================================

.. .. image:: _static/mrsimulator_cover_light.png
..     :align: center

.. only:: html

    .. cssclass:: table-bordered table-striped centered

    .. list-table::
      :widths: 25 75
      :header-rows: 0

      * - Deployment
        - .. image:: https://img.shields.io/pypi/v/mrsimulator.svg?style=flat&logo=pypi&logoColor=white
            :target: https://pypi.python.org/pypi/mrsimulator
            :alt: PyPI version

      * - Build Status
        - .. image:: https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master&logo=travis&color=white
              :target: https://travis-ci.org/DeepanshS/mrsimulator
              :alt: Build Status

          .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
              :target: https://mrsimulator.readthedocs.io/en/stable/?badge=stable
              :alt: Documentation Status

      * - License
        - .. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
              :target: https://opensource.org/licenses/BSD-3-Clause
              :alt: License

      * - Metrics
        - .. image:: https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18
              :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python
              :alt: Language grade: Python

          .. image:: https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg
              :target: https://codecov.io/gh/DeepanshS/mrsimulator

          .. image:: https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18
              :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/
              :alt: Total alerts

      * - GitHub
        - .. image:: https://img.shields.io/github/contributors/DeepanshS/mrsimulator.svg?style=flat&logo=github
              :target: https://github.com/DeepanshS/mrsimulator/graphs/contributors
              :alt: GitHub contributors

          .. image:: https://img.shields.io/github/issues-raw/deepanshs/mrsimulator
              :target: https://github.com/DeepanshS/mrsimulator/issues
              :alt: GitHub issues

.. .. image:: https://img.shields.io/pypi/dm/mrsimulator.svg?style=flat&logo=pypi
..     :target: https://img.shields.io/pypi/dm/mrsimulator
..     :alt: PyPI - Downloads

----

About
^^^^^

The `Mrsimulator` package is a C and python based library for simulating fast
solid-state nuclear magnetic resonance (NMR) spectra. The core methods are
written in C, wrapped, and made available in python for python users.

.. warning::
    The package is currently under development. We advice using with
    caution. Bug report are greatly appreciated.

----

**Features**

Fast-simulation of one-dimensional NMR spectra of uncoupled spins featuring:

- Spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}`,
  (See the list of supported isotopes),
- Arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin frequency,
- Variable angle spinning (VAS) at arbitrary angle and spin frequency,
- Static line-shape.

----

    **View the example gallery**

    .. only:: html

        .. image:: https://img.shields.io/badge/View-Example%20Gallery-Purple?s=small
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
   understanding_system
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
    ^^^^^^^^^^^^^^^^^^

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
