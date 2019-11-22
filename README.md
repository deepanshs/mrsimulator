[![Build Status](https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master)](https://travis-ci.org/DeepanshS/mrsimulator)
[![Documentation Status](https://readthedocs.org/projects/mrsimulator/badge/?version=stable)](https://mrsimulator.readthedocs.io/en/stable/?badge=stable)
[![PyPI version](https://badge.fury.io/py/mrsimulator.svg)](https://badge.fury.io/py/mrsimulator)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python)
[![codecov](https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/DeepanshS/mrsimulator)

# Mrsimulator

`Mrsimulator` is a python package with methods and tools for simulating
fast solid-state nuclear magnetic resonance (NMR) line-shapes. The
core library is written in C, wrapped and made available
in python for python users.

> :warning: The package is currently under development. We advice using with caution. Bug report are greatly appreciated.

**Features**
At present, _mrsimulator_ offers fast-simulation of one-dimensional NMR
line-shape of uncoupled spins for the following scenarios:

- Spin $I=\frac{1}{2}$, and quadrupole $I \ge \frac{1}{2}$,
  (See the list of supported isotopes),
- Arbitrary macroscopic magnetic flux density,
- Magic angle spinning (MAS) at arbitrary spin frequency,
- Variable angle spinning (VAS) at arbitrary angle and spin frequency,
- Static line-shape.

For more information, refer to the
[documentation](https://deepanshs.github.io/mrsimulator/).

## Installation

We recommend installing [anaconda](https://www.anaconda.com/distribution/)
distribution for python version 3.6 or higher. The anaconda distribution
ships with numerous packages and modules including Numpy, Scipy, and Matplotlib
which are useful for scientific datasets.

### Using pip

Pip is the probably the easiest way to install python packages.
We recommend using pip for installing Mrsimulator. PIP is a package manager
for Python packages, and is included with python version 3.4 and higher.

**Mac and Linux system**
For mac and linux systems, we provide binary distributions of mrsimulator
package for python version 3.6-3.8. To install, type the following in the
terminal,

    $ pip install mrsimulator

**Windows**
For windows, we currently do not provide binary distribution. You'll need to compile and build the mrsimulator library. Please follow the instructions below.

1.  Install dependencies using conda

        $ conda install -c conda-forge openblas fftw>=3.3.0 numpy>=1.13.3 cython>=0.29.11

2.  Install mrsimulator using pip

         $ pip install mrsimulator

### Building from source

#### Download mrsimulator package

First, clone and download the `mrsimulator` package from Github. If you prefer
`git`, type the following in the terminal,

    $ git clone git://github.com/DeepanshS/mrsimulator.git

else, [download](https://github.com/DeepanshS/mrsimulator) the package using
the browser. Once downloaded, use the terminal to navigate to the directory
containing the package (usually, the folder is named mrsimulator).

    $ cd mrsimulator

#### Installing dependencies

Next, install the dependencies of the package. We recommend using `conda` to
install the dependencies, as follows,

    $ conda install -c conda-forge openblas --file requirements.txt

#### Building and Installing mrsimulator package

Because the core of the `mrsimulator` package is written in C, you will
require a C-compiler to build and install the package. We recommend the
following C-compiler for the OS types:

- Mac OS - `clang`
- Linux - `gcc`
- Windows - `msvc`

Install the package with,

    $ python setup.py install

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/DeepanshS/mrsimulator-test/master/test_file.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://raw.githubusercontent.com/DeepanshS/mrsimulator/master/docs/_static/test_output.png)
