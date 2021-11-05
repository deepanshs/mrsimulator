.. _install:

For the users
=============

.. note::

  If you encounter an issue during instillation, see our `troubleshooting section <_troubleshooting>`__.
  If that doesn't resolve your issue, please create a bug report on our
  `Github issue tracker <https://github.com/deepanshs/mrsimulator/issues>`_.

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
    Thanks to `Sphinx-copybutton <https://sphinx-copybutton.readthedocs.io/en/latest/>`_

.. code-block:: shell

      $ python --version

For *Mac* users, Python version 3 is installed under the name ``python3``. You may replace
``python`` for ``python3`` in the above command and all subsequent python statements.

For *Windows* users, Python is not usually installed by default. See
`Python.org <https://www.python.org/downloads/windows/>`_ for a list of official Python downloads
for Windows and instillation instructions.

.. You can find more information under the Windows tab in the
.. :ref:`building_from_source` section.

.. seealso::

  If you do not have python or have an older version of python, you may visit the
  `Python downloads <https://www.python.org/downloads/>`_ on Python's official site
  and follow their installation instructions.

  If you are having issues with Python, see our
  `Python troubleshooting section <_installing_python>`__.

.. We recommend installing `Anaconda <https://www.anaconda.com/distribution/>`_
.. distribution for python version 3.6 or higher. The Anaconda distribution
.. ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
.. which are useful packages for scientific datasets. Anaconda also acts as an environment
.. manager for Python which helps keep dependencies for different projects separate.

Installing ``mrsimulator`` using pip
------------------------------------

On Google Colab Notebook
''''''''''''''''''''''''

Colaboratory is a Google research project providing a cloud-based Jupyter notebook environment.
Launch a new notebook on
`Colab <http://colab.research.google.com>`_.

Colab has an older version of numpy installed by default and it must be updated before
using mrsimulator. In a new cell, run

.. code-block:: shell

      !pip install -U numpy

and restart the runtime. Now, to install the mrsimulator package, type

.. code-block:: shell

      !pip install mrsimulator

in the first cell, and execute. All done! You may now start using the library, or
proceed to section :ref:`Introduction to Spin Systems <spin_systems_intro>` to continue
the tutorial.


.. _on_local_machine:

On Local machine (Using pip)
''''''''''''''''''''''''''''

PIP is a package manager for Python packages and is included with python version 3.4
and higher. PIP is the easiest way to install python packages.

.. tabs::

  .. tab:: Linux
    :tabid: linux

    For *Linux* users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.9. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

  .. tab:: Mac OSX
    :tabid: macosx

    For *Mac* users, we provide the binary distributions of the mrsimulator package for
    python versions 3.6-3.9. Install the package using pip as follows,

    .. code-block:: bash

        $ pip install mrsimulator

    If the above statement didn't work, you are probably using mac OS system python, in
    which case, use the following,

    .. code-block:: bash

        $ python3 -m pip install mrsimulator --user

  .. tab:: Windows
    :tabid: windows

    TODO find if we provide binary distros for windows

    .. note:: We currently do not provide binary distributions for windows. You'll need
      to compile and build the mrsimulator library from source. The following instructions
      are one-time installation only. If you are upgrading the package, see the
      :ref:`upgrading_to_a_newer_version` sub-section.

    .. include:: source_install/windows.rst

    **Install the package**.

    From within the ``Anaconda Prompt``, build and install the mrsimulator package
    using pip.

    .. code-block:: bash

      $ pip install mrsimulator

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


All done! You may now start using the library, or proceed to section
:ref:`Introduction to Spin Systems <spin_systems_intro>` to continue the tutorial.


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
`here <https://raw.githubusercontent.com/deepanshs/mrsimulator-examples/master/test_file_v0.3.py?raw=true>`_.

.. code-block:: text

    $ python test_file.py

The above statement should produce the following figure.

.. plot:: ../pyplot/test_file.py

    A test example simulation of solid-state NMR spectrum.

.. note::

    If you encounter the following error

    .. code-block:: shell

        ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject

    run

    .. code-block:: shell

        $ pip install -U numpy

    to update numpy. Running ``python test_file.py`` should now succeed.
