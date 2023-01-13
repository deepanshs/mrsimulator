
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
        - .. image:: https://img.shields.io/github/actions/workflow/status/deepanshs/mrsimulator/continuous-integration-pip.yml?logo=GitHub
            :target: https://github.com/deepanshs/mrsimulator/actions/workflows/continuous-integration-pip.yml/badge.svg
            :alt: GitHub Workflow Status

          .. image:: https://readthedocs.org/projects/mrsimulator/badge/?version=stable
            :target: https://mrsimulator.readthedocs.io/en/stable/
            :alt: Documentation Status

      * - License
        - .. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
            :target: https://opensource.org/licenses/BSD-3-Clause
            :alt: License

      * - Metrics
        - .. image:: https://codecov.io/gh/deepanshs/mrsimulator/branch/master/graph/badge.svg
            :target: https://codecov.io/gh/deepanshs/mrsimulator

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

.. only:: html

  **About**

.. only:: not html

  About
  '''''

**mrsimulator** is an open-source Python package for fast computation/analysis of nuclear
magnetic resonance (NMR) spectra in fluid and solid phases.


----

.. only:: html

    .. raw:: html

        <h3>See our example galleries</h3>
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
- It is a fast and versatile multi-dimensional solid-state NMR spectra simulator, including MAS
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- It includes simulation of weakly coupled nuclei experiencing J and dipolar couplings.
- It is fully documented with a stable and simple API and is easily incorporated into
  Python scripts and web apps.
- It is compatible with modern Python packages, such as Scikit-learn, Keras, etc.
- Packages using **mrsimulator** -

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
    ax[0].plot(sim.methods[0].simulation)
    ax[0].set_title("Static")
    ax[1].plot(sim.methods[1].simulation)
    ax[1].set_title("MAS")
    plt.tight_layout()
    plt.show()

----

.. note::

  Throughout the web version of this documentation, you can copy code blocks into your clipboard by
  hovering over the top right corner of each gray code block and clicking the copy-to-clipboard
  icon. This is useful for copying code examples into your Python scripts and Jupyter notebooks.


**Features**

The **mrsimulator** package offers the following

- **Fast simulation** of one and two-dimensional solid-state NMR spectra.

- Simulation of **coupled and uncoupled spin system**
    - for spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei
    - at arbitrary macroscopic magnetic flux density
    - at arbitrary rotor angles
    - at arbitrary spinning frequency

- A library of pre-built **NMR methods**,
    - 1D Bloch decay spectrum
    - 1D Bloch decay central transition spectrum
    - 2D Multi-Quantum Variable Angle Spinning (MQ-VAS)
    - 2D Satellite-Transition Variable Angle Spinning (ST-VAS)
    - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT)
    - 2D Magic-Angle Flipping (MAF)
    - 2D Dynamic-Angle Spinning (DAS)
    - Custom user-defined methods (Method)

- **Models** for tensor parameter distribution in amorphous materials.
    - Czjzek
    - Extended Czjzek
    - Custom user-defined models

----

Introduction
------------

.. toctree::
    :maxdepth: 2
    :caption: Introduction

    installation/installation
    introduction/getting_started
    introduction/isotopomers_example
    introduction/fitting_example

User Guide
----------

.. toctree::
    :maxdepth: 3
    :caption: User Documentation

    user_guide/spin_system/spin_system
    user_guide/spin_system_distributions/spin_system_distributions
    user_guide/methods_library/methods_library
    user_guide/method/method
    user_guide/simulator/simulator
    user_guide/signal_processor/signal_processor
    user_guide/io/mrsim_IO

Examples
--------

.. toctree::
    :maxdepth: 2
    :caption: Galleries

    examples/index
    fitting/index
    signal_processor/index

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
    api_py/py-signal-processor
    api_py/py-model
    api_py/py-fitting

.. api_c/c_api


Project details
---------------

.. toctree::
    :maxdepth: 1
    :caption: Project details

    changelog
    credits/contributors
    credits/license
    credits/acknowledgment


Reporting Bugs
--------------

Submit bug reports or feature requests on the
`Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_.

Discussions are welcome on the `Github discussion <https://github.com/deepanshs/mrsimulator/discussions>`_ page.


How to cite
-----------

Please refer to `mrsimulator Github page <https://github.com/deepanshs/mrsimulator>`_ for details.

.. only:: html

    Index
    -----

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
