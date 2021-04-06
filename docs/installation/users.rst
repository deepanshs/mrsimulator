.. _install:

For the users
=============

Strict Requirements
-------------------

``mrsimulator`` has the following strict requirements:

- `Python <https://www.python.org>`_ 3.6 or later
- `Numpy <https://numpy.org>`_ 1.17 or later

See :ref:`requirements` for a full list of requirements.

Make sure you have the required version of python by typing the following in the terminal,

.. tip::
    You may also click the copy-button located at the top-right corner of the code cell
    area in the HTML docs, to copy the code lines without the prompts and then paste it
    as usual.
    Thanks to `Sphinx-copybutton <https://sphinx-copybutton.readthedocs.io/en/latest/>`_)

.. code-block:: shell

      $ python --version

For *Mac* users, python version 3 is installed under the name *python3*. You may replace
*python* for *python3* in the above command and all subsequent python statements.

For *Windows* users, we recommend the `Anaconda <https://www.anaconda.com/products/individual/>`_
or `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ distribution of
python>3.6. Anaconda distribution for python comes with popular python packages that
are frequently used in scientific computing.
Miniconda is a minimal installer for conda. It is a smaller version of Anaconda that
includes conda, Python, and the packages they depend on, along with other useful
packages such as pip.

.. You can find more information under the Windows tab in the
.. :ref:`building_from_source` section.

.. seealso::

  If you do not have python or have an older version of python, you may visit the
  `Python downloads <https://www.python.org/downloads/>`_ or
  `Anaconda <https://www.anaconda.com/products/individual/>`_ websites and follow their
  instructions on how to install python.

.. We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
.. distribution for python version 3.6 or higher. The anaconda distribution
.. ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
.. which are useful packages for scientific datasets.

Installing ``mrsimulator`` using pip
------------------------------------

On Google Colab Notebook
''''''''''''''''''''''''

Colaboratory is a Google research project. It is a Jupyter notebook environment that
runs entirely in the cloud. Launch a new notebook on
`Colab <http://colab.research.google.com>`_. To install the mrsimulator package, type

.. code-block:: shell

      !pip install mrsimulator

in the first cell, and execute. All done! You may now proceed to the next section and
start using the library.



.. _on_local_machine:

On Local machine (Using pip)
''''''''''''''''''''''''''''

PIP is a package manager for Python packages and is included with python version 3.4
and higher. PIP is the easiest way to install python packages.

.. tabs::

  .. tab:: Linux
    :tabid: linux

    For *Linux* users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.8. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

  .. tab:: Mac OSX
    :tabid: macosx

    For *Mac* users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.8. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

    If the above statement didn't work, you are probably using mac OS system python, in
    which case, use the following,

    .. code-block:: bash

        $ python3 -m pip install mrsimulator --user

  .. tab:: Windows
    :tabid: windows

    .. note:: We currently do not provide binary distributions for windows. You'll need
      to compile and build the mrsimulator library from source. The following instructions
      are one-time installation only. If you are upgrading the package, see the
      :ref:`upgrading_to_a_newer_version` sub-section.

    .. include:: windows.rst

    **Install the package**.

    From within the ``Anaconda Prompt``, build and install the mrsimulator package
    using pip.

    .. code-block:: bash

      $ pip install mrsimulator

If you get a ``PermissionError``, it usually means that you do not have the required
administrative access to install new packages to your Python installation. In this
case, you may consider adding the ``--user`` option, at the end of the statement, to
install the package into your home directory. You can read more about how to do this in
the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.

.. _upgrading_to_a_newer_version:

Upgrading to a newer version
""""""""""""""""""""""""""""

If you are upgrading to a newer version of ``mrsimulator``, you have all the prerequisites
installed on your system. In this case, type the following in the terminal/Prompt

.. code-block:: bash

    $ pip install mrsimulator -U


All done! You may now proceed to the next section and start using the library.


.. _building_from_source:

Building from the source
------------------------

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
`PyPI <https://pypi.org/project/mrsimulator/#files>`_.


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

    **OpenBLAS and FFTW libraries**

    On Linux, the package manager for your distribution is usually the easiest route to
    ensure you have the prerequisites to building the mrsimulator library. To build from
    source, you will need the OpenBLAS and FFTW development headers for your Linux
    distribution. Type the following command in the terminal, based on your Linux
    distribution.

    *For (Debian/Ubuntu):*

    .. code-block:: bash

      $ sudo apt-get install libopenblas-dev libfftw3-dev

    *For (Fedora/RHEL):*

    .. code-block:: bash

      $ sudo yum install openblas-devel fftw-devel

    **Install a C/C++ compiler**

    The C-compiler comes with your Linux distribution. No further action is
    required.

  .. tab:: Mac OSX
    :tabid: macosx_source

    **OpenBLAS/Accelerate and FFTW libraries**

    You will require the ``brew`` package manager to install the development headers for the
    OpenBLAS (if applicable) and FFTW libraries. Read more on installing brew from
    `homebrew <https://brew.sh>`_.

    *Step-1* Install the FFTW library using the `homebrew <https://brew.sh>`_ formulae.

    .. code-block:: bash

      $ brew install fftw

    *Step-2* By default, the mrsimulator package links to the openblas library for BLAS
    operations. Mac users may opt to choose the in-build Apple's Accelerate library. If you
    opt for Apple's Accelerate library, skip to `Step-3`. If you wish to link the mrsimulator
    package to the OpenBLAS library, type the following in the terminal,

    .. code-block:: bash

      $ brew install openblas

    *Step-3* If you choose to link the mrsimulator package to the OpenBLAS library, skip
    to the next section, :ref:`building_and_installing`.

    *(a)* You will need to install the BLAS development header for Apple's Accelerate
    library. The easiest way is to install the Xcode Command Line Tools. Note, this is a
    one-time installation. If you have previously installed the Xcode Command Line Tools,
    you may skip this sub-step. Type the following in the terminal,

    .. code-block:: bash

      $ xcode-select --install

    *(b)* The next step is to let the mrsimulator setup know your preference.
    Open the ``settings.py`` file, located at the root level of the mrsimulator source
    code folder, in a text editor. You should see

    .. code-block:: python

      # -*- coding: utf-8 -*-
      # BLAS library
      use_openblas = True
      # mac-os only
      use_accelerate = False

    To link the mrsimulator package to the Apple's Accelerate library, change the
    fields to

    .. code-block:: python

      # -*- coding: utf-8 -*-
      # BLAS library
      use_openblas = False
      # mac-os only
      use_accelerate = True

    **Install a C/C++ compiler**

    The C-compiler installs with the Xcode Command Line Tools. No further action is
    required.

  .. tab:: Windows
    :tabid: windows_source

    .. include:: windows.rst


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
`here <https://raw.githubusercontent.com/DeepanshS/mrsimulator-examples/master/test_file_v0.3.py?raw=true>`_.

.. code-block:: text

    $ python test_file.py

The above statement should produce the following figure.

.. plot:: ../pyplot/test_file.py

    A test example simulation of solid-state NMR spectrum.
