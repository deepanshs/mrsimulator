

.. _shielding_tensor_api:

============
Installation
============

Installing requirements
+++++++++++++++++++++++

``mrsimulator`` requires `NumPy>=1.13.3 <http://www.numpy.org>`_
and intel `mkl <https://pypi.org/project/mkl/>`_ and
`mkl_include <https://pypi.org/project/mkl-include>`_ C routines to build
and install the package. You can install these dependent libraries using

.. code-block:: text

    pip install -r requirements.txt


A c compiler is required to successful compile and build the ``mrsimulator``
package.

On linux, you may install gcc c-compiler.

.. code-block:: text

    sudo apt install gcc

Installing mrsimulator
++++++++++++++++++++++

To install the ``mrsimulator`` package, type the following
in the terminal.

.. code-block:: text

    pip install git+https://github.com/DeepanshS/mrsimulator.git@master


Test and verify the build
+++++++++++++++++++++++++

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
