

.. _shielding_tensor_api:

================================
Installing `mrsimulator` package
================================

We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
distribution for python version 3.6 or higher. The anaconda distribution
ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
which are useful packages for scientific datasets.

Downloading mrsimulator package
*******************************

First, clone and download the ``mrsimulator`` package from Github. If you use
``git``, type the following in the terminal,

.. code-block:: bash

    $ git clone git://github.com/DeepanshS/mrsimulator.git

else, `download <https://github.com/DeepanshS/mrsimulator>`_ the package using
your browser. Once downloaded, use the terminal to change the current directory
location to the downloaded folder (usually, the downloader folder is named
mrsimulator).

.. code-block:: bash

    $ cd mrsimulator


Installing dependencies
***********************

Next, install the dependencies of the package. We recommend using ``conda`` to
install the dependencies, as follows,

.. code-block:: bash

    $ conda install -c conda-forge openblas --file requirements.txt



Building and Installing mrsimulator package
*******************************************

Since the core of the ``mrsimulator`` package is written in C, you will require
a C compiler to build and install the package.

.. code-block:: shell

    $ python setup.py install

.. pip install git+https://github.com/DeepanshS/mrsimulator.git@master



Check your build
****************

If the installation is successful, you should be able to run the following test
file in your terminal. Download the test file `here <link>`_.

.. code-block:: text

    $ python mrsimulator_quick_test.py

This will display the following message on the screen

.. code-block:: text

    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 0.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.

    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 0.0 ppm
    Shielding anisotropy = 13.89 ppm
    Shielding asymmetry = 0.25
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 1000.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.

    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 0.0 ppm
    Shielding anisotropy = 13.89 ppm
    Shielding asymmetry = 0.25

and the corresponding plot shown below.

.. image:: /_static/test_output.png
