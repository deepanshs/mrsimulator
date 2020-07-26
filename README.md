# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://badge.fury.io/py/mrsimulator.svg)](https://badge.fury.io/py/mrsimulator) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mrsimulator)                                                                                                                                                                                                                                                                                                                                                                                                              |
| Build Status | [![Travis (.org)](https://img.shields.io/travis/deepanshs/mrsimulator?logo=travis)](https://travis-ci.org/github/DeepanshS/mrsimulator) [![GitHub Workflow Status](<https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI%20(pip)?logo=GitHub>)](https://github.com/DeepanshS/mrsimulator/actions) [![Read the Docs](https://img.shields.io/readthedocs/mrsimulator)](https://mrsimulator.readthedocs.io/en/latest/)                                                                                                                                                       |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| Metrics      | [![PyPI - Downloads](https://img.shields.io/pypi/dm/mrsimulator.svg)](https://img.shields.io/pypi/dm/mrsimulator) [![Total alerts](https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python) [![codecov](https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/DeepanshS/mrsimulator) |

`mrsimulator` is a python package for computing fast real-time solid-state
magnetic resonance (NMR) line-shapes/spectrum. The library is optimized to compute bulk
solid-state line-shapes, enabling the simulation of both crystalline and amorphous-like
materials. The core of the `mrsimulator` library is written in C, wrapped and made
available in python.

> :warning: The package is currently under development. We advice using with caution. Bug report are greatly appreciated.

## Features

At present, the `mrsimulator` package offers the following

- **Real-time simulation** of one-dimensional solid-state NMR line-shapes. See our
  [benchmark results](https://mrsimulator.readthedocs.io/en/latest/benchmark.html).

- **Uncoupled spin-system**

  - for spin I=1/2, and quadrupole I>1/2 nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- The library includes the following **NMR methods**,

  - 1D Bloch decay spectrum, and
  - 1D Bloch decay central transition spectrum.

## Goals for the near future

Our current objectives for the future are the following

- Include line-shape simulation of coupled spin-systems for

  - spin I=1/2, and quadrupole I>1/2 nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- Expand the library of NMR methods. We expect to include the following methods

  - 2D Multi-Quantum Magic Angle Spinning (MQ-MAS),
  - 2D Dynamic Angle Spinning (DAS),
  - 2D Magic Angle Flipping (MAF), and
  - 2D isotropic to anisotropic sideband correlation spectrum (PASS).

For more information, refer to the
[documentation](https://mrsimulator.readthedocs.io/en/latest/).

> **View our example gallery**
>
> [![](https://img.shields.io/badge/View-Example%20Gallery-Purple?s=small)](https://mrsimulator.readthedocs.io/en/latest/auto_examples/index.html)

## Installation

    $ pip install mrsimulator

Please read our [installation document](https://mrsimulator.readthedocs.io/en/latest/installation.html) for details.

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/DeepanshS/mrsimulator-examples/master/test_file_v0.3.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://mrsimulator.readthedocs.io/en/master/_images/test_file.png)
