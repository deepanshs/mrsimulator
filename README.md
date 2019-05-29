
# mrsimulator

The package ``mrsimulator`` contains key functionality and tools needed for
simulating one-dimensional nuclear magnetic resonance (NMR) lineshapes using
Python. The package is currently under development. Version 0.1.0 supports
simulation of single spin \(I=1/2\) nuclei, static,
magic angle spinning (MAS), and variable angle spinning (VSA) lineshapes.
For more information, read the documentation [here](https://deepanshs.github.io/mrsimulator/).

## Installation

### Installing requirements

The `mrsimulator` package requires the [fftw3](https://anaconda.org/eumetsat/fftw3)
C routines to function. Download and install the fftw3 by typing the following
in the terminal

    conda install -c eumetsat fftw3

Additionally, ``mrsimulator`` requires the following packages.

- [NumPy>=1.13.3](http://www.numpy.org) for array manipulation.
- [astropy>=3.0](https://www.astropy.org) for units library.
- mkl, mkl-include for lineshape calculation.

<!-- You may install these package using pip as follows

    pip install numpy scipy astropy mkl mkl-include -->

For figures and visualization, we use

- [matplotlib>=3.0.2](https://matplotlib.org)

For the web-face interface, we use

- [plotly>=3.6](https://plot.ly/python/)
- dash>=0.40
- dash_daq>=0.1

<!-- To install these package with

    pip install matplotlib plotly dash dash_daq -->

### Installing mrsimulator

To install the ``mrsimulator`` package, type the following
in the terminal.

    pip install git+https://github.com/DeepanshS/mrsimulator.git@master

## Test

If the installation is successful, you should be able to run the following
in the terminal.

    python -c "import mrsimulator; mrsimulator.run_test()"

This will display the following message on the screen

    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 0.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.

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

    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    isotropic chemical shift = 0.0 ppm
    chemical shift anisotropy = 13.89 ppm
    chemical shift asymmetry = 0.25

and a corresponding plot shown below.

![alt text](https://raw.githubusercontent.com/DeepanshS/mrsimulator/gh-pages/_static/test_output.png)
