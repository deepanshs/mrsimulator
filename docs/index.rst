
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
        - .. image:: https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI%20(pip)?logo=GitHub
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

          .. image:: https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge
            :target: https://www.codefactor.io/repository/github/deepanshs/mrsimulator
            :alt: CodeFactor

      * - GitHub
        - .. image:: https://img.shields.io/github/contributors/DeepanshS/mrsimulator.svg?logo=github
            :target: https://github.com/DeepanshS/mrsimulator/graphs/contributors
            :alt: GitHub contributors

          .. image:: https://img.shields.io/github/issues/deepanshs/mrsimulator?logo=github
            :target: https://github.com/DeepanshS/mrsimulator/issues
            :alt: GitHub issues

      * - Citation
        - .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3978779.svg
            :target: https://doi.org/10.5281/zenodo.3978779

.. .. image:: https://img.shields.io/github/commits-since/deepanshs/mrsimulator/v0.2.1?logo=github
..   :alt: GitHub commits since tagged version

.. .. image:: https://img.shields.io/pypi/dm/mrsimulator.svg?style=flat&logo=pypi
..     :target: https://img.shields.io/pypi/dm/mrsimulator
..     :alt: PyPI - Downloads

**About**

``mrsimulator`` is an open-source python package for fast simulation and analysis of
multi-dimensional solid-state magnetic resonance (NMR) spectra of crystalline materials,
bio macro-molecules, and amorphous materials. Simulate the NMR spectrum of macro-molecules
or amorphous in just a few seconds.

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
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding, paramagnetic shift)
  and quadrupolar coupling interactions.
- Includes simulation of weakly coupled nuclei experiencing J and dipolar couplings.
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

- The library of included **NMR methods** follow,
    - 1D Bloch decay spectrum, and
    - 1D Bloch decay central transition spectrum.
    - 2D Multi-quantum Variable Angle Spinning (MQ-VAS),
    - 2D Satellite-transition Variable Angle Spinning (ST-VAS), and
    - 2D Dynamic Angle Spinning (DAS),
    - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT), and
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
    api_py/py-fitting
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
`Github issue tracker <https://github.com/DeepanshS/mrsimulator/issues>`_. Reports
are also welcomed  by directly contacting `Deepansh Srivastava <mailto:srivastava.89@osu.edu>`_.

Discussions are welcome on `Github discussion <https://github.com/DeepanshS/mrsimulator/discussions>`_


.. only:: html

    Index
    -----

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
