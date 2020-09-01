
########################################
Welcome to the Mrsimulator documentation
########################################

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

          .. image:: https://img.shields.io/pypi/pyversions/mrsimulator
            :alt: PyPI - Python Version

      * - Build Status
        - .. image:: https://img.shields.io/travis/deepanshs/mrsimulator?logo=travis
            :target: https://travis-ci.org/DeepanshS/mrsimulator
            :alt: Travis (.org)

          .. image:: https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI%20(pip)?logo=GitHub
            :target: https://github.com/DeepanshS/mrsimulator/actions
            :alt: GitHub Workflow Status

          .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
            :target: https://mrsimulator.readthedocs.io/en/stable/
            :alt: Documentation Status

      * - License
        - .. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
            :target: https://opensource.org/licenses/BSD-3-Clause
            :alt: License

      * - Metrics
        - .. image:: https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm
            :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python
            :alt: Language grade: Python

          .. image:: https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg
            :target: https://codecov.io/gh/DeepanshS/mrsimulator

          .. image:: https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm
            :target: https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/
            :alt: Total alerts

      * - GitHub
        - .. image:: https://img.shields.io/github/contributors/DeepanshS/mrsimulator.svg?logo=github
            :target: https://github.com/DeepanshS/mrsimulator/graphs/contributors
            :alt: GitHub contributors

          .. image:: https://img.shields.io/github/issues/deepanshs/mrsimulator?logo=github
            :target: https://github.com/DeepanshS/mrsimulator/issues
            :alt: GitHub issues

      * - Citation
        - .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3978780.svg
            :target: https://doi.org/10.5281/zenodo.3978780

.. .. image:: https://img.shields.io/github/commits-since/deepanshs/mrsimulator/v0.2.1?logo=github
..   :alt: GitHub commits since tagged version

.. .. image:: https://img.shields.io/pypi/dm/mrsimulator.svg?style=flat&logo=pypi
..     :target: https://img.shields.io/pypi/dm/mrsimulator
..     :alt: PyPI - Downloads

**About**

``mrsimulator`` is an open-source python package for fast computation/analysis of solid-state
magnetic resonance (NMR) spectra of both crystalline and amorphous materials. The core of
the ``mrsimulator`` library is written in C, wrapped, and made available in python.

----

**Why use mrsimulator?**

- It is open-source and free.
- It is a fast and versatile solid-state NMR simulator of one-dimensional static, MAS,
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- Future release will include simulations of weakly coupled nuclei experiencing J and dipolar
  couplings, and multi-dimensional NMR spectra.
- It is fully documented with a stable and simple API and is easily incorporated into your
  python scripts and web apps.
- It is compatible with modern python packages, such as scikit-learn, Keras, etc.
- Packages using mrsimulator -

  - `mrinversion <https://mrinversion.readthedocs.io/en/latest/>`_

----

.. only:: html

    .. raw:: html

        <br>

    **View our example gallery**

    .. image:: https://img.shields.io/badge/View-Example%20Gallery-Purple?s=small
        :target: auto_examples/index.html

----

**Features**

The ``mrsimulator`` package currently offers the following

- **Fast simulation** of one-dimensional solid-state NMR spectra. See our :ref:`benchmark` results.

- Simulation of **uncoupled spin system**
    - for spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei,
    - at arbitrary macroscopic magnetic flux density,
    - at arbitrary rotor angles, and
    - at arbitrary spinning frequency.

- The library includes the following **NMR methods**,
    - 1D Bloch decay spectrum, and
    - 1D Bloch decay central transition spectrum.

- **Models** for tensor parameter distribution in amorphous materials.
    - Czjzek
    - Extendend Czjzek

----

**Goals for the near future**

Our current objectives are the following

- Include spectral simulation of coupled spin systems for
    - spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei,
    - at arbitrary macroscopic magnetic flux density,
    - at arbitrary rotor angles, and
    - at arbitrary spinning frequency.

- Expand the library of NMR methods. We expect to include the following methods
    - 2D Multi-quantum Magic Angle Spinning (MQ-MAS),
    - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT).
    - 2D Dynamic Angle Spinning (DAS), and
    - 2D Magic Angle Flipping (MAF).


.. Contribution
.. ------------

.. ``Mrsimulator`` is a open source NMR simulation package. We are a small team
.. working on developing the package for the NMR community. Any contribution and
.. suggestion is greatly appreciated.


.. warning::
    The package is currently under development. We advice using with
    caution. Bug report are greatly appreciated.

----

Getting Started
---------------

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    installation
    requirements
    introduction
    getting_started
    getting_started-objects
    configuring_simulator
    mrsim_IO

Signal Processing (``mrsimulator.SignalProcessor``)
---------------------------------------------------

.. toctree::
    :maxdepth: 2
    :caption: Signal Processing

    signal_processing

Models
------

.. toctree::
    :maxdepth: 2
    :caption: Models

    model/czjzek
    model/extended_czjzek


Examples and Benchmarks
-----------------------

.. toctree::
    :maxdepth: 2
    :caption: Examples and Benchmarks

    auto_examples/index
    benchmark

Theory
------

.. toctree::
    :maxdepth: 2
    :caption: Theory

    theory/components
    theory/models

API and references
------------------

.. toctree::
    :maxdepth: 2
    :caption: API and references

    api_py/py-simulator
    api_py/py-signal-processing
    api_py/py-model
    api_c/c_api


Project details
---------------

.. toctree::
    :maxdepth: 1
    :caption: Project details

    changelog
    credits/contributors
    credits/license
    credits/acknowledgment

.. .. only:: html

.. 	.. toctree::
.. 		:maxdepth: 2
.. 		:caption: Table of Contents:

.. 		about
.. 		installation
.. 		requirements
.. 		understanding-spin_systems
.. 		getting_started
.. 		using_mrsimulator_objects
.. 		load_sample
.. 		configuring_simulator
.. 		benchmark
.. 		auto_examples/index
.. 		theory/components
.. 		api_py/py_api
.. 		api_c/c_api


.. understanding_system
..    objects
..    spectrum_object
..    theory/wigner_rotations
..    examples


.. only:: html

    Index
    -----

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
