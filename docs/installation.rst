

.. _shielding_tensor_api:

================================
Installing `Mrsimulator` package
================================

We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
distribution for python version 3.6 or higher. The anaconda distribution
ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
which are useful packages for scientific datasets.

Using pip
---------

Pip is the probably the easiest way to install python packages. PIP is a
package manager for Python packages, and is included with python version 3.4
and higher.

**Mac and Linux system**

For mac and linux systems, we provide binary distributions of mrsimulator
package for python version 3.6 and 3.7. To install, type the following in the
terminal,

.. code-block:: bash

    $ pip install mrsimulator

**Windows or python version>=3.8**

For windows, or python version>=3.8, you need to compile and build the
mrsimulator library. Follow the instructions below.

1) Install dependencies using conda
    .. code-block:: bash

        $ conda install -c conda-forge openblas fftw>=3.3.0 numpy>=1.13.3 cython>=0.29.11

2) Install mrsimulator using pip
    .. code-block:: bash

        $ pip install mrsimulator


From source
-----------

**Download mrsimulator package**

`Clone and download <https://github.com/DeepanshS/mrsimulator>`_ the mrsimulator
package from Github. If you use git, type the following in the terminal,

.. code-block:: bash

    $ git clone git://github.com/DeepanshS/mrsimulator.git

Once downloaded, use the terminal to navigate to the directory
containing the package (usually, the folder is named mrsimulator).

.. code-block:: bash

    $ cd mrsimulator


**Installing dependencies**

Next, install the dependencies of the package. We recommend using conda to
install the dependencies, as follows,

.. code-block:: bash

    $ conda install -c conda-forge openblas --file requirements.txt



**Building and Installing mrsimulator package**

Because the core of the mrsimulator package is written in C, you will
require a C-compiler to build and install the package. We recommend the
following C-compiler for the OS types:

- Mac OS - ``clang``
- Linux - ``gcc``
- Windows - ``msvc``

Install the package with,

.. code-block:: shell

    $ python setup.py install

.. pip install git+https://github.com/DeepanshS/mrsimulator.git@master


Check your build
----------------

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file
`here <https://raw.github.com/DeepanshS/mrsimulator-test/master/test_file.py?raw=true>`_.

.. code-block:: text

    $ python test_file.py

This should produce the following figure.

.. figure:: _static/test_output.*
    :figclass: figure-polaroid
