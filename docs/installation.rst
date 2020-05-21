

.. _shielding_tensor_api:

================================
Installing `Mrsimulator` package
================================

.. We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
.. distribution for python version 3.6 or higher. The anaconda distribution
.. ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
.. which are useful packages for scientific datasets.

Using pip
---------

PIP is a package manager for Python packages, and is included with python version 3.4
and higher. PIP is the simplest way to install python packages.

**Mac and Linux system**

For mac and linux systems, we provide binary distributions of mrsimulator
package for python version 3.6-3.8. To install, type the following in the
terminal,

.. code-block:: bash

    $ pip install mrsimulator

**Windows**

We currently do not provide binary distribution. You'll need to
compile and build the mrsimulator library. Follow the instructions below.

1. Install dependencies using conda

.. code-block:: bash

    $ conda install -c conda-forge openblas fftw>=3.3.0 numpy>=1.13.3 cython>=0.29.11

2. Build package using pip (you will require a c-compiler, such as,
`Microsoft Visual C++ <https://visualstudio.microsoft.com/downloads/>`_).

.. code-block:: bash

    $ pip install mrsimulator


Building from source
--------------------
.. The binary distribution of the mrsimulator package includes the above two libraries.

Mrsimulator utilizes the BLAS and FFTW libraries for computation. To leverage the best
performance from the mrsimulator package, we recommend installing the BLAS and FFTW
libraries, which may be optimized for your system. In the following, we
list a few recommendations on how to install the BLAS, FFTW, and mrsimulator libraries.


Download the mrsimulator package
''''''''''''''''''''''''''''''''

`Clone and download <https://github.com/DeepanshS/mrsimulator>`_ the mrsimulator
package from Github. If you prefer git, type the following in the terminal,

.. code-block:: bash

    $ git clone git://github.com/DeepanshS/mrsimulator.git

Once downloaded, use the terminal to navigate to the directory
containing the package (usually, the folder is named mrsimulator).

.. code-block:: bash

    $ cd mrsimulator


Installation
''''''''''''

MacOS users
***********

**Installing dependencies**

**Step-1** By default, the mrsimulator package links to the openblas library for BLAS
operations. Mac users may opt to choose the in-build apple's accelerate library. If you
opt for apple's accelerate library, skip to step-2. If you wish to link the mrsimulator
package to the openblas library, follow

.. code-block:: bash

    $ brew install openblas

**Step-2** Install the FFTW library using the `homebrew <https://brew.sh>`_ formulae,
and the remaining dependencies using pip,

.. code-block:: bash

    $ brew install fftw
    $ pip install -r requirements.txt

.. $ conda install -c conda-forge openblas --file requirements.txt

**Building and installing the mrsimulator package**

Because the core of the mrsimulator package is written in C, you will
require a C-compiler to build and install the package.

**Step-3** If you choose to link the
mrsimulator package to openblas library, skip to step-4.
Open the ``setting.py`` file, which is located at the root level of the mrsimulator
folder. You should see,

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

**Step-4** Install the package.

.. code-block:: bash

    $ python setup.py install

.. pip install git+https://github.com/DeepanshS/mrsimulator.git@master


Linux(Ubuntu) users
*******************

**Installing dependencies**

**Step-1** For Ubuntu users, openblas and FFTW libraries may already be installed. If
not, install the libraries with

.. code-block:: bash

    $ sudo apt-get install libopenblas-dev libfftw3-dev

**Step-2** Install the remaining dependencies using pip.

.. code-block:: bash

    $ pip install -r requirements.txt

**Building and installing the mrsimulator package**

**Step-3** Install the package.

.. code-block:: bash

    $ python setup.py install

Linux(CentOS) users
*******************

**Installing dependencies**

**Step-1** Install the openblas and FFTW libraries.

.. code-block:: bash

    $ yum install openblas-devel fftw-devel

**Step-2** Install the remaining dependencies using pip.

.. code-block:: bash

    $ pip install -r requirements.txt

**Building and installing the mrsimulator package**

**Step-3** Install the package.

.. code-block:: bash

    $ python setup.py install

.. We recommend the
.. following C-compiler for the OS types:

.. - Mac OS - ``clang``
.. - Linux - ``gcc``
.. - Windows - ``msvc``

Check your build
----------------

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file
`here <https://raw.github.com/DeepanshS/mrsimulator-test/master/test_file_v0.3.py?raw=true>`_.

.. code-block:: text

    $ python test_file.py

This should produce the following figure.

.. figure:: _static/test_output.*
    :figclass: figure-polaroid
