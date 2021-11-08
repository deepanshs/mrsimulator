
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
            :target: #
            :alt: PyPI - Python Version

      * - Build Status
        - .. image:: https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI?logo=GitHub
            :target: https://github.com/deepanshs/mrsimulator/actions
            :alt: GitHub Workflow Status

          .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
            :target: https://mrsimulator.readthedocs.io/en/stable/
            :alt: Documentation Status

      * - License
        - .. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
            :target: https://opensource.org/licenses/BSD-3-Clause
            :alt: License

      * - Metrics
        - .. image:: https://img.shields.io/lgtm/grade/python/g/deepanshs/mrsimulator.svg?logo=lgtm
            :target: https://lgtm.com/projects/g/deepanshs/mrsimulator/context:python
            :alt: Language grade: Python

          .. image:: https://codecov.io/gh/deepanshs/mrsimulator/branch/master/graph/badge.svg
            :target: https://codecov.io/gh/deepanshs/mrsimulator

          .. image:: https://img.shields.io/lgtm/alerts/g/deepanshs/mrsimulator.svg?logo=lgtm
            :target: https://lgtm.com/projects/g/deepanshs/mrsimulator/alerts/
            :alt: Total alerts

          .. image:: https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge
            :target: https://www.codefactor.io/repository/github/deepanshs/mrsimulator
            :alt: CodeFactor

      * - Social
        - .. image:: https://img.shields.io/github/contributors/deepanshs/mrsimulator?style=social&logo=github
            :target: https://github.com/deepanshs/mrsimulator/graphs/contributors
            :alt: GitHub contributors

          .. image:: https://img.shields.io/github/issues/deepanshs/mrsimulator?style=social&logo=github
            :target: https://github.com/deepanshs/mrsimulator/issues
            :alt: GitHub issues

          .. image:: https://img.shields.io/github/stars/deepanshs/mrsimulator?style=social
            :target: https://github.com/deepanshs/mrsimulator/stargazers
            :alt: GitHub stars

          .. image:: https://img.shields.io/github/forks/deepanshs/mrsimulator?style=social
            :target: https://github.com/deepanshs/mrsimulator/network/members
            :alt: GitHub forks

.. .. image:: https://img.shields.io/github/commits-since/deepanshs/mrsimulator/v0.2.1?logo=github
..   :alt: GitHub commits since tagged version

.. .. image:: https://img.shields.io/pypi/dm/mrsimulator.svg?style=flat&logo=pypi
..     :target: https://img.shields.io/pypi/dm/mrsimulator
..     :alt: PyPI - Downloads

**About**

``mrsimulator`` is an open-source python package for fast simulation and analysis of
multi-dimensional solid-state magnetic resonance (NMR) spectra of crystalline and
amorphous materials.

----

.. only:: html

    .. raw:: html

        <h3>See our example gallery</h3>
        <div class='sim-gallery'>
          <div>
            <a href="examples/index.html">Gallery
              <p></p>
              Simulation Examples
            </a>
          </div>
          <div>
            <a href="fitting/index.html">Gallery
              <p></p>
              <a href="fitting/index.html">Fitting Examples
            </a>
          </div>
        </div>

----

**Why use mrsimulator?**

- It is open-source and free.
- It is a fast and versatile multi-dimensional solid-state NMR spectra simulator including, MAS
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- It includes simulation of weakly coupled nuclei experiencing J and dipolar couplings.
- It is fully documented with a stable and simple API and is easily incorporated into your
  python scripts and web apps.
- It is compatible with modern python packages, such as scikit-learn, Keras, etc.
- Packages using mrsimulator -

  - `mrinversion <https://mrinversion.readthedocs.io/en/stable/>`_

----


**Features**

The ``mrsimulator`` package offers the following

- **Fast simulation** of one/two-dimensional solid-state NMR spectra. See our :ref:`benchmark` results.

- Simulation of **coupled and uncoupled spin system**
    - for spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei,
    - at arbitrary macroscopic magnetic flux density,
    - at arbitrary rotor angles, and
    - at arbitrary spinning frequency.

- A library of **NMR methods**,
    - 1D Bloch decay spectrum,
    - 1D Bloch decay central transition spectrum,
    - 2D Multi-quantum Variable Angle Spinning (MQ-VAS),
    - 2D Satellite-transition Variable Angle Spinning (ST-VAS),
    - 2D Dynamic Angle Spinning (DAS),
    - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT),
    - 2D Magic Angle Flipping (MAF).

- **Models** for tensor parameter distribution in amorphous materials.
    - Czjzek
    - Extended Czjzek

----

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

    installation/installation
    introduction
    getting_started
    getting_started-objects
    getting_started_ethanol
    configuring_simulator
    mrsim_IO
.. designing_methods

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

    examples/index
    fitting/index
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
    api_py/py-utilities
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
.. 		examples/index
.. 		theory/components
.. 		api_py/py_api
.. 		api_c/c_api


.. understanding_system
..    objects
..    spectrum_object
..    theory/wigner_rotations
..    examples

Reporting Bugs
--------------

The preferred location for submitting feature requests and bug reports is the
`Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_. Reports
are also welcomed  by directly contacting `Deepansh Srivastava <mailto:srivastava.89@osu.edu>`_.

Discussions are welcome on `Github discussion <https://github.com/deepanshs/mrsimulator/discussions>`_


How to cite
-----------

If you use mrsimulator in your publication, please consider citing the following.

- Deepansh J. Srivastava, Maxwell Venetos, Philip J. Grandinetti, Shyam Dwaraknath, & Alexis McCarthy. (2021, May 26). mrsimulator: v0.6.0 (Version v0.6.0). Zenodo. http://doi.org/10.5281/zenodo.4814638

- Srivastava DJ, Vosegaard T, Massiot D, Grandinetti PJ (2020) Core Scientific Dataset Model: A lightweight and portable model and file format for multi-dimensional scientific data. PLOS ONE 15(1): e0225953. https://doi.org/10.1371/journal.pone.0225953

*Additionally, if you use lmfit for least-squares fitting, consider citing the lmfit package.*

- Matt Newville; Renee Otten; Andrew Nelson; Antonino Ingargiola; Till Stensitzki; Dan Allan; Austin Fox; Faustin Carter; Michał; Dima Pustakhod; lneuhaus; Sebastian Weigand; Ray Osborn; Glenn; Christoph Deil; Mark; Allan L. R. Hansen; Gustavo Pasquevich; Leon Foks; Nicholas Zobrist; Oliver Frost; Alexandre Beelen; Stuermer; kwertyops; Anthony Polloreno; Shane Caldwell; Anthony Almarza; Arun Persaud; Ben Gamari; Benjamin F. Maier. (2021, February 7). lmfit/lmfit-py 1.0.2 (Version 1.0.2). Zenodo. http://doi.org/10.5281/zenodo.4516651

.. only:: html

    Index
    -----

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
