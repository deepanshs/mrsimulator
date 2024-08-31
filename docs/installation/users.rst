.. _install:

For the users
=============

.. note::

  If you encounter an issue during installation, see our
  `troubleshooting section <instillation_troubleshooting>`_.
  If that doesn't resolve your issue, please create a bug report on our
  `Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_.\


Strict Requirements
-------------------

**MRSimulator** has the following strict requirements:

- `Python <https://www.python.org>`_ |PY_VERSION| or later
- `Numpy <https://numpy.org>`_ 1.17 or later

See :ref:`requirements` for a full list of requirements.

Make sure you have the required version of Python by typing the following in the terminal,

.. code-block:: shell

      $ python --version

For *MacOS* users, Python version 3 is installed under the name *python3*. You may
replace *python* for *python3* in the above command and all subsequent Python statements.

For *Windows* users, Python is not usually installed by default. See
`Python.org <https://www.python.org/downloads/windows/>`_ for a list of official Python
downloads and Windows installation instructions.

.. seealso::

  If you do not have Python or have an older version of Python, you may visit the
  `Python downloads <https://www.python.org/downloads/>`_ or
  `Anaconda <https://www.anaconda.com/products/individual/>`_ websites and follow their
  instructions on installing Python.

Installing mrsimulator
----------------------

.. only:: html

  .. tabs::

    .. tab:: Google Colab Notebook

      .. include:: colab.rst

    .. tab:: Local machine (Using pip)

      .. include:: pip.rst

    .. tab:: From source

      .. include:: source.rst

.. only:: not html

  Google Colab Notebook
  '''''''''''''''''''''
  .. include:: colab.rst

  Local machine (Using pip)
  '''''''''''''''''''''''''
  .. include:: pip.rst

  From source
  '''''''''''
  .. include:: source.rst


Updating mrsimulator
--------------------


If you are upgrading to a newer version of mrsimulator, you should have all the prerequisites
already installed on your system. In this case, type the following in the terminal/Prompt

.. code-block:: bash

    $ pip install mrsimulator -U


All done! You may now start using the library or proceed to
:ref:`getting_started` to continue the tutorial.

Testing your build
------------------

.. note::
  For Windows users using anaconda Python 3.9 and higher, you need to set the following
  environment variable in the ``Anaconda Prompt`` before running mrsimulator scripts.

  .. code-block:: bash

      $ set CONDA_DLL_SEARCH_MODIFICATION_ENABLE='1'

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file
`here <https://raw.githubusercontent.com/deepanshs/mrsimulator-examples/master/test_file_v0.3.py?raw=true>`_
or copy and paste the following code into a Python file and run the code.

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
    ax[0].plot(sim.methods[0].simulation.real)
    ax[0].set_title("Static")
    ax[1].plot(sim.methods[1].simulation.real)
    ax[1].set_title("MAS")
    plt.tight_layout()
    plt.show()

.. note::

    If you encounter the following error

    .. code-block:: shell

        ValueError: numpy.ndarray size changed, may indicate binary incompatibility.
        Expected 88 from C header, got 80 from PyObject

    update numpy by running

    .. code-block:: shell

        $ pip install -U numpy
