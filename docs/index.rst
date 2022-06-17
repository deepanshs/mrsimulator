
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
              Simulation
            </a>
          </div>
          <div>
            <a href="fitting/index.html">Gallery
              <p></p>
              <a href="fitting/index.html">
              Fitting
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

**A brief example**

.. skip: next

.. plot::
    :caption: Simulation of static and MAS solid-state NMR spectra

    from mrsimulator import Simulator, SpinSystem, Site
    from mrsimulator.method.lib import BlochDecaySpectrum
    import matplotlib.pyplot as plt

    # Make Site and SpinSystem objects
    H_site = Site(isotope="1H", shielding_symmetric={"zeta": 13.89, "eta": 0.25})
    spin_system = SpinSystem(sites=[H_site])

    # Make static and MAS one-pulse acquire Method objects
    static = BlochDecaySpectrum(channels=["1H"])
    mas = BlochDecaySpectrum(channels=["1H"], rotor_frequency=1000)  # in Hz

    # Setup and run the Simulation object
    sim = Simulator(spin_systems=[spin_system], methods=[static, mas])
    sim.run()

    # Plot the spectra
    fig, ax = plt.subplots(1, 2, figsize=(6, 3), subplot_kw={"projection": "csdm"})
    ax[0].plot(sim.methods[0].simulation.real, color="black", linewidth=1)
    ax[0].set_title("Static")
    ax[1].plot(sim.methods[1].simulation.real, color="black", linewidth=1)
    ax[1].set_title("MAS")
    plt.tight_layout()
    plt.show()


----


**Features**

The ``mrsimulator`` package offers the following

- **Fast simulation** of one and two-dimensional solid-state NMR spectra.

- Simulation of **coupled and uncoupled spin system**
    - for spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei
    - at arbitrary macroscopic magnetic flux density
    - at arbitrary rotor angles
    - at arbitrary spinning frequency

- A library of **NMR methods**,
    - 1D Bloch decay spectrum
    - 1D Bloch decay central transition spectrum
    - 2D Multi-quantum Variable Angle Spinning (MQ-VAS)
    - 2D Satellite-transition Variable Angle Spinning (ST-VAS)
    - 2D Dynamic Angle Spinning (DAS)
    - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT)
    - 2D Magic Angle Flipping (MAF)
    - Custom user-defined 1D and 2D methods (Method)

- **Models** for tensor parameter distribution in amorphous materials.
    - Czjzek
    - Extended Czjzek
    - Custom user-defined models

----

.. Contribution
.. ------------

.. ``Mrsimulator`` is a open source NMR simulation package. We are a small team
.. working on developing the package for the NMR community. Any contribution and
.. suggestion is greatly appreciated.

----

Introduction
------------

.. toctree::
    :maxdepth: 2
    :caption: Introduction

    installation/installation
    introduction/getting_started
    introduction/ethanol_example

User Guide
----------

.. toctree::
    :maxdepth: 2
    :caption: User Documentation

    user_guide/spin_system/spin_system
    user_guide/spin_system_distributions/spin_system_distributions
    user_guide/methods_library/methods_library
    user_guide/method/method_advanced_user
    user_guide/method/query_objects
    user_guide/simulator/simulator
    user_guide/signal_processing/signal_processing
    user_guide/io/mrsim_IO

Examples
--------

.. toctree::
    :maxdepth: 2
    :caption: Galleries

    examples/index
    fitting/index
    signal_processing/index

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

Submit bug reports or feature requests on the
`Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_.

Discussions are welcome on the `Github discussion <https://github.com/deepanshs/mrsimulator/discussions>`_ page.


How to cite
-----------

If you use mrsimulator in your publication, please consider citing the following.

- Deepansh J. Srivastava, Matthew Giammar, Maxwell C. Venetos, Shyam Dwaraknath, Philip J. Grandinetti, & Alexis McCarthy. (2021). mrsimulator: v0.6.1. Zenodo. https://doi.org/10.5281/zenodo.5559730

- Srivastava DJ, Vosegaard T, Massiot D, Grandinetti PJ (2020) Core Scientific Dataset Model: A lightweight and portable model and file format for multi-dimensional scientific data. PLOS ONE 15(1): e0225953. https://doi.org/10.1371/journal.pone.0225953

*Additionally, if you use lmfit for least-squares fitting, consider citing the lmfit package.*

- Matt Newville; Renee Otten; Andrew Nelson; Antonino Ingargiola; Till Stensitzki; Dan Allan; Austin Fox; Faustin Carter; Michał; Dima Pustakhod; lneuhaus; Sebastian Weigand; Ray Osborn; Glenn; Christoph Deil; Mark; Allan L. R. Hansen; Gustavo Pasquevich; Leon Foks; Nicholas Zobrist; Oliver Frost; Alexandre Beelen; Stuermer; kwertyops; Anthony Polloreno; Shane Caldwell; Anthony Almarza; Arun Persaud; Ben Gamari; Benjamin F. Maier. (2021, February 7). lmfit/lmfit-py 1.0.2 (Version 1.0.2). Zenodo. https://doi.org/10.5281/zenodo.4516651

.. only:: html

    Index
    -----

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
