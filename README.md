
[![Build Status](https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master)](https://travis-ci.org/DeepanshS/mrsimulator)
[![Documentation Status](https://readthedocs.org/projects/mrsimulator/badge/?version=stable)](https://mrsimulator.readthedocs.io/en/stable/?badge=stable)

# mrsimulator

``mrsimulator`` is a library package with methods and tools for fast
simulation of solid-state nuclear magnetic resonance (NMR) line-shapes.
The library contains routines written in C which are wrapped and made
available in python.

The package is currently under development. At present, `mrsimulator` features
simulation of one-dimensional NMR line-shape of uncoupled spin
I=1/2 isotopes for the following scenarios --

- At arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin rate,
- Variable angle spinning (VAS) at arbitrary angle and spin rates,
- Static line-shape.

For more information, refer to the
[documentation](https://deepanshs.github.io/mrsimulator/).

## Installation

We recommend installing [anaconda](https://www.anaconda.com/distribution/)
distribution for python version 3.6 or higher. The anaconda distribution
ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
which are useful packages for scientific datasets.

### Download mrsimulator package

First, clone and download the ``mrsimulator`` package from Github. If you prefer
``git``, type the following in the terminal,

    $ git clone git://github.com/DeepanshS/mrsimulator.git

else, [download](https://github.com/DeepanshS/mrsimulator) the package using
the browser. Once downloaded, use the terminal to navigate to the directory
containing the package (usually, the folder is named mrsimulator).

    $ cd mrsimulator

### Installing dependencies

Next, install the dependencies of the package. We recommend using ``conda`` to
install the dependencies, as follows,

    $ conda install -c conda-forge openblas --file requirements.txt

### Building and Installing mrsimulator package

Because the core of the ``mrsimulator`` package is written in C, you will
require a C-compiler to build and install the package. We recommend the
following C-compiler for the OS types:

- Mac OS - ``clang``
- Linux - ``gcc``
- Windows - ``msvc``

Install the package with,

    $ python setup.py install

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/DeepanshS/mrsimulator-test/master/test_file.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://raw.githubusercontent.com/DeepanshS/mrsimulator/master/docs/_static/test_output.png)
