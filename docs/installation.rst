

.. _install:

================================
Installing `Mrsimulator` package
================================

The `Mrsimulator` package requires python version 3.6 or higher. Make sure you have
python>=3.6 installed by typing the following in the terminal (for Mac and Linux
users only),

.. code-block:: shell

      $ python --version

For `Mac` users, python version 3 is installed under the name `python3`. You may replace
`python` for `python3` in the above command.

For `Windows` users, follow the instruction under the Windows tab in the section
:ref:`on_local_machine`.

.. seealso::

  If you do not have python or have an older version of python, you may visit
  `Python <https://www.python.org/downloads/>`_ or
  `anaconda <https://www.anaconda.com/products/individual/>`_ and follow their
  instructions to install python.

.. We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
.. distribution for python version 3.6 or higher. The anaconda distribution
.. ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
.. which are useful packages for scientific datasets.

.. _on_local_machine:

On Local machine
----------------

.. PIP is a package manager for Python packages, and is included with python version 3.4
.. and higher. PIP is the simplest way to install python packages.

.. tabs::

  .. tab:: Linux

    For Linux users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.8. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

  .. tab:: Mac OSX

    For `Mac` users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.8. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

    If the above statement didn't work, you are probably using mac os system python, in
    which case, use the following,

    .. code-block:: bash

        $ python3 -m pip install mrsimulator --user

  .. tab:: Windows

    .. note:: We currently do not provide binary distributions for windows. You'll need
        to compile and build the mrsimulator library. Follow the instructions below.

    .. include:: install-docs/windows.rst

    **Step-4**. Build and install the package.

    From within the ``Anaconda Prompt``, build and install the mrsimulator package
    using pip.

    .. code-block:: bash

      $ pip install mrsimulator

----

On Google Colab Notebook
------------------------

Colaboratory is a Google research project. It is a Jupyter notebook environment that
runs entirely in the cloud. Launch a new notebook on
`Colab <http://colab.research.google.com>`_. To install the mrsimulator package, type

.. code-block:: shell

      !pip install mrsimulator

in the first cell, and execute. All done! You may now start using the library.

.. _building_from_source:

Building from the source
------------------------

The Mrsimulator package utilizes the BLAS and FFTW libraries for computation. To
leverage the best performance from the mrsimulator package, we recommend installing the
BLAS and FFTW libraries, which may be optimized for your system. In the following, we
list a few recommendations on how to install the BLAS, FFTW, and mrsimulator libraries.

**Step-A**: Clone and download the mrsimulator package from
`Github <https://github.com/DeepanshS/mrsimulator>`_.

**Step-B**: Install OS-specific dependencies.

.. tabs::

  .. tab:: Linux (Ubuntu)

    For Ubuntu users, openblas and FFTW libraries may already be installed with the OS.
    If not, install the libraries with

    .. code-block:: bash

      $ sudo apt-get install libopenblas-dev libfftw3-dev

  .. tab:: Linux (CentOS)

    Install the openblas and FFTW libraries.

    .. code-block:: bash

      $ yum install openblas-devel fftw-devel

  .. tab:: Mac OSX

    **Step-1** By default, the mrsimulator package links to the openblas library for BLAS
    operations. Mac users may opt to choose the in-build apple's accelerate library. If you
    opt for apple's accelerate library, skip to Step-2. If you wish to link the mrsimulator
    package to the openblas library, install openblas using `homebrew <https://brew.sh>`_
    formulae as follows,

    .. code-block:: bash

      $ brew install openblas

    **Step-2** Install the FFTW library using the `homebrew <https://brew.sh>`_ formulae.

    .. code-block:: bash

      $ brew install fftw

    **Step-3** If you choose to link the mrsimulator package to openblas library, skip
    to Step-C. Open the ``settings.py`` file, located at the root level of the
    mrsimulator folder, in a text editor. You should see,

    .. code-block:: python

      # -*- coding: utf-8 -*-
      # BLAS library
      use_openblas = True
      # mac-os only
      use_accelerate = False

    To link the mrsimulator package to the in-build apple's accelerate library, change the
    fields to

    .. code-block:: python

      # -*- coding: utf-8 -*-
      # BLAS library
      use_openblas = False
      # mac-os only
      use_accelerate = True

  .. tab:: Windows

    .. include:: install-docs/windows.rst

**Step-C**: Use the terminal/Prompt to navigate into the directory containing the
package (usually, the folder is named mrsimulator),

.. code-block:: bash

    $ cd mrsimulator

and install the package using pip,

.. code-block:: bash

    $ pip install .



.. We recommend the
.. following C-compiler for the OS types:

.. - Mac OS - ``clang``
.. - Linux - ``gcc``
.. - Windows - ``msvc`` (https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019)


----

Setup for developers and contributors
-------------------------------------

A GitHub account is required for developers and contributors. Make sure you have
git installed on your system.

**Step-A** (Optional) Create a virtual environment. It is a good practice to create
separate virtual python environments for packages when in developer mode.
The following is an example of a Conda environment.

.. code-block:: bash

    $ conda create -n mrsimulator-dev python=3.7
    $ conda activate mrsimulator-dev

**Step-B** Clone the mrsimulator repository using git and navigate into the package
folder.

.. code-block:: bash

    $ git clone git://github.com/DeepanshS/mrsimulator.git
    $ cd mrsimulator

**Step-C** Follow the instruction under Step-B from :ref:`building_from_source`
section. For developers and contributors using mac OSX, please run the setup
by binding to the openblas libraries.

**Step-D** Build and install the package in the development mode using pip.

.. code-block:: bash

    $ pip install -e .

**Step-E**: Install the requirements for developers using pip.

.. code-block:: bash

    $ pip install -r requirements-dev.txt

Note for the developers and contributors
''''''''''''''''''''''''''''''''''''''''

**Running tests**: For unit tests, we use the pytest module. At the root directory
of the mrsimulator package folder, type

.. code-block:: bash

    $ pytest

which will run a series of tests.

**Building docs**: We use the sphinx python documentation generator for building docs.
Navigate to the ``docs`` folder within the mrsimulator package folder, and type,

.. code-block:: bash

    $ make html

The above command will build the documentation and store the build at
``mrsimulator/docs/_build/html``. Double click the `index.html` file within this
folder to view the offline documentation.

**Submitting pull requests** Make sure all the test pass and the documentation build
is successful before creating a pull request.

.. We recommend the
.. following C-compiler for the OS types:
.. - Mac OS - ``clang``
.. - Linux - ``gcc``
.. - Windows - ``msvc`` (https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019)



Check your build
----------------

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file
`here <https://raw.githubusercontent.com/DeepanshS/mrsimulator-examples/master/test_file_v0.3.py?raw=true>`_.

.. code-block:: text

    $ python test_file.py

The above statement should produce the following figure.

.. figure:: _static/test_output.*
    :figclass: figure
