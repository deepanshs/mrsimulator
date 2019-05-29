

.. _shielding_tensor_api:

============
Installation
============

Installing requirements
+++++++++++++++++++++++

The `mrsimulator` package requires `fftw3 <https://anaconda.org/eumetsat/fftw3>`_
C routines. Download and install the fftw3 routines by typing
the following in the terminal

.. code-block:: python

    conda install -c eumetsat fftw3

In addition, `mrsimulator` also requires `NumPy>=1.13.3 <http://www.numpy.org>`_
and intel `mkl <https://pypi.org/project/mkl/>`_ and
`mkl_include <https://pypi.org/project/mkl-include/>`_ C routines to build and
install the `mrsimulator` package. Download and install these C routines using

.. code-block:: python

    pip install "numpy>=1.13.1" mkl mkl-include

Some additional package dependencies are

 - `astropy>=3.0 <https://www.astropy.org>`_ for the units library,
 - `matplotlib>=3.0.2 <https://matplotlib.org>`_ for figures and visualization,

and,

 - `plotly>=3.6 <https://plot.ly/python/>`_
 - `dash>=0.40 <https://pypi.org/project/dash/>`_
 - `dash_daq>=0.1 <https://pypi.org/project/dash-daq/>`_

for the web-face interface.


Installing mrsimulator
++++++++++++++++++++++

To install the ``mrsimulator`` package, type the following
in the terminal.

.. code-block:: python

    pip install git+https://github.com/DeepanshS/mrsimulator.git@master


Test
++++

If the installation is successful, you should be able to run the following
in the terminal.

.. code-block:: python

    python -c "import mrsimulator; mrsimulator.run_test()"

This will display the following message on the screen

.. code-block:: python

    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 0.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 0.0 ppm
    chemical shift anisotropy = 13.89 ppm
    chemical shift asymmetry = 0.25
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 1000.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 0.0 ppm
    chemical shift anisotropy = 13.89 ppm
    chemical shift asymmetry = 0.25

and the corresponding plot shown below.

.. image:: /_static/test_output.png
