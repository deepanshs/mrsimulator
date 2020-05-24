# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://badge.fury.io/py/mrsimulator.svg)](https://badge.fury.io/py/mrsimulator)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| Build Status | [![Build Status](https://travis-ci.org/DeepanshS/mrsimulator.svg?branch=master)](https://travis-ci.org/DeepanshS/mrsimulator) [![Documentation Status](https://readthedocs.org/projects/mrsimulator/badge/?version=master)](https://mrsimulator.readthedocs.io/en/master/?badge=master)                                                                                                                                                                                                                                                                                                                              |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| Metrics      | [![PyPI - Downloads](https://img.shields.io/pypi/dm/mrsimulator.svg)](https://img.shields.io/pypi/dm/mrsimulator)[![Total alerts](https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python) [![codecov](https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/DeepanshS/mrsimulator) |

`Mrsimulator` is a python package for computing fast real-time solid-state nuclear
magnetic resonance (NMR) line-shapes/spectrum. The library is optimized to compute bulk
solid-state line-shapes, enabling the simulation of both crystalline and amorphous-like
materials. The core of the `Mrsimulator` library is written in C, wrapped and made
available in python.

> :warning: The package is currently under development. We advice using with caution. Bug report are greatly appreciated.

## Features

At present, the `Mrsimulator` package offers a fast real-time simulation of
one-dimensional solid-state NMR line-shape of

- Uncoupled spin-system for

  - spin I=1/2, and quadrupole I>1/2` nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- The library includes the following NMR methods,

  - 1D Bloch decay spectrum, and
  - 1D Bloch decay central transition spectrum.

## Goals for the near future

Our current objectives for the future are the following

- Include line-shape simulation of coupled spin-systems for

  - spin :math:`I=\frac{1}{2}`, and quadrupole :math:`I \ge \frac{1}{2}` nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- Expand the library of NMR methods. We expect to include the following methods

  - 2D Multi-Quantum Magic Angle Spinning (MQ-MAS),
  - 2D Dynamic Angle Spinning (DAS),
  - 2D Magic Angle Flipping (MAF), and
  - 2D isotropic to anisotropic sideband correlation spectrum (PASS).

For more information, refer to the
[documentation](https://mrsimulator.readthedocs.io/en/stable/).

> **View our example gallery**
>
> [![](https://img.shields.io/badge/View-Example%20Gallery-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/auto_examples/index.html)

## Installation

PIP is a package manager for Python packages, and is included with python version 3.4
and higher. PIP is the simplest way to install python packages.

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

Mrsimulator utilizes the BLAS and FFTW libraries for computation. To leverage the best
performance from the mrsimulator package, we recommend installing the BLAS and FFTW
libraries, which may be optimized for your system. In the following, we
list a few recommendations on how to install the BLAS, FFTW, and mrsimulator libraries.

#### Download mrsimulator package

Clone and [download](https://github.com/DeepanshS/mrsimulator) the mrsimulator package
from Github. If you prefer `git`, type the following in the terminal,

    $ git clone git://github.com/DeepanshS/mrsimulator.git

Once downloaded, use the terminal to navigate to the directory
containing the package (usually, the folder is named mrsimulator).

    $ cd mrsimulator

#### Installation

##### MacOS users

###### Installing dependencies

**Step-1** By default, the mrsimulator package links to the openblas library for BLAS
operations. Mac users may opt to choose the in-build apple's accelerate library. If you
opt for apple's accelerate library, skip to step-2. If you wish to link the mrsimulator
package to the openblas library, follow

    $ brew install openblas

**Step-2** Install the FFTW library using the [homebrew](https://brew.sh) formulae,
and the remaining dependencies using pip,

    $ brew install fftw
    $ pip install -r requirements.txt

###### Building and Installing mrsimulator package

Because the core of the mrsimulator package is written in C, you will
require a C-compiler to build and install the package.

**Step-3** If you choose to link the
mrsimulator package to openblas library, skip to step-4.
Open the `setting.py` file, which is located at the root level of the mrsimulator
folder. You should see,

    # -*- coding: utf-8 -*-
    # BLAS library
    use_openblas = True
    # mac-os only
    use_accelerate = False

To link the mrsimulator package to the in-build apple's accelerate library, change the
fields to

    # -*- coding: utf-8 -*-
    # BLAS library
    use_openblas = False
    # mac-os only
    use_accelerate = True

**Step-4** Install the package.

    $ python setup.py install

##### Linux(Ubuntu) users

###### Installing dependencies

**Step-1** For Ubuntu users, openblas and FFTW libraries may already be installed. If
not, install the libraries with

    $ sudo apt-get install libopenblas-dev libfftw3-dev

**Step-2** Install the remaining dependencies using pip.

    $ pip install -r requirements.txt

###### Building and installing the mrsimulator package

**Step-3** Install the package.

    $ python setup.py install

##### Linux(CentOS) users

###### Installing dependencies

**Step-1** Install the openblas and FFTW libraries.

    $ yum install openblas-devel fftw-devel

**Step-2** Install the remaining dependencies using pip.

    $ pip install -r requirements.txt

###### Building and installing the mrsimulator package

**Step-3** Install the package.

    $ python setup.py install

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/DeepanshS/mrsimulator-test/master/test_file_v0.3.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://raw.githubusercontent.com/DeepanshS/mrsimulator/master/docs/_static/test_output.png)
