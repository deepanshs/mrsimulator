
**Prerequisites**

You will need a C-compiler suite and the development headers for the BLAS and FFTW libraries,
along with development headers from Python and Numpy, to build the
**MRSimulator** library from source.
The **MRSimulator** package utilizes the BLAS and FFTW routines for numerical computation.
To leverage the best performance, we recommend installing the BLAS and FFTW libraries which are
optimized and tuned for your system. In the following, we list recommendations on installing the
C-compiler (if applicable), BLAS, FFTW, and building the **MRSimulator** libraries.

**Obtaining the Source Packages**

The latest stable source package for **MRSimulator** is available on
`PyPI <https://pypi.org/project/mrsimulator/#files>`_ and
`Github release <https://github.com/deepanshs/mrsimulator/releases>`_. Download and
extract the *.tar.gz* file.

**OS-dependent prerequisites**

.. note::
 Installing OS-dependent prerequisites is a one-time process. If upgrading to a newer version of
 mrsimulator, skip to the next section.

.. only:: html

  .. tabs::

    .. tab:: Linux

      .. include:: source_install/linux.rst

    .. tab:: Mac OSX

      .. include:: source_install/macosx.rst

    .. tab:: Windows

      .. include:: source_install/windows.rst


.. only:: not html

  **Linux**

  .. include:: source_install/linux.rst

  **Mac OSX**

  .. include:: source_install/macosx.rst

  **Windows**

  .. include:: source_install/windows.rst


**Building and Installing**

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
