.. _install:

For the users
=============

.. note::

   If you encounter an issue during installation, see our
   `troubleshooting section <instillation_troubleshooting>`_.

   If that doesn't resolve your issue, please create a bug report on our
   `Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_.

Strict Requirements
-------------------

``mrsimulator`` has the following strict requirements:

- `Python <https://www.python.org>`_ 3.6 or later
- `Numpy <https://numpy.org>`_ 1.17 or later

See :ref:`requirements` for a full list of requirements.

Make sure you have the required version of python by typing the following in the terminal,

.. code-block:: shell

      $ python --version

For *Mac* users, python version 3 is installed under the name *python3*. You may replace
*python* for *python3* in the above command and all subsequent python statements.

For *Windows* users, Python is not usually installed by default. See
`Python.org <https://www.python.org/downloads/windows/>`_ for a list of official Python downloads
and Windows installation instructions.

.. seealso::

  If you do not have python or have an older version of python, you may visit the
  `Python downloads <https://www.python.org/downloads/>`_ or
  `Anaconda <https://www.anaconda.com/products/individual/>`_ websites and follow their
  instructions on how to install python.

Installing ``mrsimulator``
--------------------------

On Google Colab Notebook
''''''''''''''''''''''''

Colaboratory is a Google research project. It is a Jupyter notebook environment that
runs entirely in the cloud. Launch a new notebook on
`Colab <http://colab.research.google.com>`_. We recommend going through the *Welcome to Colab!*
tutorial if you are new to Notebooks.

By default, Colaboratory has an older version of ``numpy`` installed which first needs to be
updated. In a new cell, run

.. code-block:: shell

    !pip install -U numpy

and press the *Restart Runtime* button

To install the ``mrsimulator`` package, type

.. code-block:: shell

    !pip install mrsimulator

in a new cell, and execute. All done! You may now start using the library, or
proceed to :ref:`getting_started` to continue the tutorial.


.. _on_local_machine:

On Local machine (Using pip)
''''''''''''''''''''''''''''

PIP is a package manager for Python packages and is included with python version 3.4
and higher. PIP is the easiest way to install python packages. Install the package
using pip as follows,

.. code-block:: bash

    $ pip install mrsimulator

For *Mac* users, if the above statement didn't work, you are probably using mac OS
system python, in which case, use the following,

.. code-block:: bash

    $ python3 -m pip install mrsimulator --user

For windows users using anaconda python 3.8 and higher, you need to set the following
environment variable in the ``Anaconda Prompt`` before running mrsimulator scripts.

.. code-block:: bash

    $ set CONDA_DLL_SEARCH_MODIFICATION_ENABLE='1'

If you get a ``PermissionError``, it usually means that you do not have the required
administrative access to install new packages to your Python installation. In this
case, you may consider adding the ``--user`` option at the end of the statement to
install the package into your home directory. You can read more about how to do this in
the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

.. _upgrading_to_a_newer_version:

Upgrading to a newer version
""""""""""""""""""""""""""""

If you are upgrading to a newer version of ``mrsimulator``, you have all the prerequisites
installed on your system. In this case, type the following in the terminal/Prompt

.. code-block:: bash

    $ pip install mrsimulator -U


All done! You may now start using the library, or proceed to
:ref:`getting_started` to continue the tutorial.


.. _building_from_source:

Building ``mrsimulator`` from the source
----------------------------------------

Prerequisites
'''''''''''''

You will need a C-compiler suite and the development headers for the BLAS and FFTW
libraries, along with development headers from Python and Numpy, to build the
``mrsimulator`` library from source.
The mrsimulator package utilizes the BLAS and FFTW routines for numerical computation.
To leverage the best performance, we recommend installing the BLAS and FFTW libraries,
which are optimized and tuned for your system. In the following,
we list recommendations on how to install the c-compiler (if applicable), BLAS, FFTW,
and building the mrsimulator libraries.

Obtaining the Source Packages
"""""""""""""""""""""""""""""

Stable packages
***************

The latest stable source package for ``mrsimulator`` is available on
`PyPI <https://pypi.org/project/mrsimulator/#files>`_ and
`Github  release <https://github.com/deepanshs/mrsimulator/releases>`_. Download and
extract the *.tar.gz* file.


.. _os_dependent_prerequisite:

OS-dependent prerequisites
""""""""""""""""""""""""""

.. note::
    Installing OS-dependent prerequisites is a one-time process. If you are
    upgrading to a newer version of mrsimulator, skip to :ref:`building_and_installing`
    section.

.. tabs::

  .. tab:: Linux
    :tabid: linus_source

    .. include:: source_install/linux.rst

  .. tab:: Mac OSX
    :tabid: macosx_source

    .. include:: source_install/macosx.rst

  .. tab:: Windows
    :tabid: windows_source

    .. include:: source_install/windows.rst


.. _building_and_installing:

Building and Installing
"""""""""""""""""""""""

Use the terminal/Prompt to navigate into the directory containing the
package (usually, the folder is named mrsimulator),

.. code-block:: bash

    $ cd mrsimulator

From within the source code folder, type the following in the terminal to install the
library.

.. code-block:: bash

    $ pip install .

If you get an error that you don't have the permission to install the package into
the default ``site-packages`` directory, you may try installing with the ``--user``
options as,

.. code-block:: bash

    $ pip install . --user


Test your build
---------------

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file
`here <https://raw.githubusercontent.com/deepanshs/mrsimulator-examples/master/test_file_v0.3.py?raw=true>`_
or copy and paste the following code into a python file and run the code.

.. skip: next

.. plot::
    :caption: An example simulating solid-state NMR spectrums of static and MAS experiments

    from mrsimulator import Simulator, SpinSystem, Site
    from mrsimulator.methods import BlochDecaySpectrum
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

.. note::

    If you encounter the following error

    .. code-block:: shell

        ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject

    update numpy by running

    .. code-block:: shell

        $ pip install -U numpy
