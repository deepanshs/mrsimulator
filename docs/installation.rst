

.. _shielding_tensor_api:

================================
Installing `mrsimulator` package
================================

We recommend installing `anaconda <https://www.anaconda.com/distribution/>`_
distribution for python version 3.6 or higher. The anaconda distribution
ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
which are useful packages for scientific datasets. In addition,
conda provides `mkl optimized <https://docs.anaconda.com/mkl-optimizations/>`_
for numerical libraries such as Numpy, Scipy.

.. If you have opted for the ``nomkl``, we suggest you create a new conda
.. environment before proceeding. You can read about creating new conda
.. environment `here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_.




Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

Clone and download the `mrsimulator` package,

.. code-block:: shell

    $ git clone git://github.com/DeepanshS/mrsimulator.git
    $ cd mrsimulator


.. and install the dependencies using

.. .. code-block:: shell

..     $ cd mrsimulator
..     $ conda install --file requirements.txt


Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

In Anaconda versions 2.5 and later, intel MKL is freely available by default.
To build and link ``mrsimulator`` with the intel-mkl libraries follow,

.. code-block:: shell

    $ conda install mkl mkl-include --file requirements.txt

If you, however, wish to opt out of MKL and instead use
`openBLAS <http://www.openblas.net/>`_, execute the following lines.

.. code-block:: shell

    $ conda install nomkl openblas --file requirements.txt

A c compiler is required to successful compile and build the ``mrsimulator``
package.

.. On linux, you can get the gcc compiler.

.. .. code-block:: text

..     $ sudo apt install gcc


To install the ``mrsimulator`` package, type the following
in the terminal.

.. code-block:: text

    python setup.py install

.. pip install git+https://github.com/DeepanshS/mrsimulator.git@master


Test
++++

If the installation is successful, you should be able to run the following
in the terminal.

.. code-block:: text

    python -c "import mrsimulator; mrsimulator.run_test()"

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
