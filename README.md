
[![Build Status](https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master)](https://travis-ci.org/DeepanshS/mrsimulator)
[![Documentation Status](https://readthedocs.org/projects/mrsimulator/badge/?version=stable)](https://mrsimulator.readthedocs.io/en/stable/?badge=stable)

# mrsimulator

The package ``mrsimulator`` contains key functionality and tools needed for
simulating one-dimensional nuclear magnetic resonance (NMR) lineshapes using
Python. The package is currently under development. Version 0.1.0 supports
simulation of single spin \(I=1/2\) nuclei, static,
magic angle spinning (MAS), and variable angle spinning (VSA) lineshapes.
For more information, read the documentation [here](https://deepanshs.github.io/mrsimulator/).

## Installation

### Installing requirements

``mrsimulator`` requires [NumPy>=1.13.3](http://www.numpy.org)
and intel [mkl](https://pypi.org/project/mkl/) and
[mkl_include](https://pypi.org/project/mkl-include) C routines to build and
install the package. You can install these dependent libraries using

    pip install -r requirements.txt

A c compiler is required to successful compile and build the ``mrsimulator``
package.

On linux

    sudo apt install gcc

### Installing mrsimulator

To install the ``mrsimulator`` package, type the following
in the terminal.

    pip install git+https://github.com/DeepanshS/mrsimulator.git@master

## Test and verify the build

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

![alt text](https://raw.githubusercontent.com/DeepanshS/mrsimulator/gh-pages/_static/test_output.png)
