# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://img.shields.io/pypi/v/mrsimulator.svg?style=flat&logo=pypi&logoColor=white)](https://pypi.python.org/pypi/mrsimulator) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mrsimulator)                                                                                                                                                                                                                                                                                                                                                                                                          |
| Build Status | [![GitHub Workflow Status](<https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI%20(pip)?logo=GitHub>)](https://github.com/DeepanshS/mrsimulator/actions) [![Read the Docs](https://img.shields.io/readthedocs/mrsimulator)](https://mrsimulator.readthedocs.io/en/stable/)                                                                                                                                                                                                                                                                                                                                         |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Metrics      | [![Total alerts](https://img.shields.io/lgtm/alerts/g/DeepanshS/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/DeepanshS/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/DeepanshS/mrsimulator/context:python) [![codecov](https://codecov.io/gh/DeepanshS/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/DeepanshS/mrsimulator) [![CodeFactor](https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge)](https://www.codefactor.io/repository/github/deepanshs/mrsimulator) |
| Citation     | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3978780.svg)](https://doi.org/10.5281/zenodo.3978779)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

`mrsimulator` is an open-source python package for fast computation/analysis of solid-state
magnetic resonance (NMR) spectra of both crystalline and amorphous materials. The core
of the `mrsimulator` library is written in C, wrapped, and made available in python.

> :warning: The package is currently under development. We advice using with caution. Bug report are greatly appreciated.

## Why use mrsimulator?

- It is open-source and free.
- It is a fast and versatile solid-state NMR simulator of one-dimensional static, MAS,
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- Include simulations of weakly coupled nuclei experiencing J and dipolar couplings.
- It is fully documented with a stable and simple API and is easily incorporated into your
  python scripts and web apps.
- It is compatible with modern python packages, such as scikit-learn, Keras, etc.
- Packages using mrsimulator -
  - [mrinversion](https://mrinversion.readthedocs.io/en/stable/)

> **View our example gallery**
>
> [![Simulation](https://img.shields.io/badge/View-Simulation%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/examples/index.html)
>
> [![Fitting](https://img.shields.io/badge/View-Fitting%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/fitting/index.html)

### Features

The `mrsimulator` package currently offers the following

- **Fast simulation** of one-dimensional solid-state NMR spectra. See our
  [benchmark results](https://mrsimulator.readthedocs.io/en/stable/benchmark.html).

- Simulation of **coupled and uncoupled spin system**

  - for spin I=1/2, and quadrupole I>1/2 nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- The library includes the following **NMR methods**,

  - 1D Bloch decay spectrum,
  - 1D Bloch decay central transition spectrum,
  - 2D Multi-Quantum Variable Angle Spinning (MQ-VAS),
  - 2D Satellite-transition Variable Angle Spinning (MQ-VAS), and
  - 2D Dynamic Angle Spinning (DAS),
  - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT), and
  - 2D Magic Angle Flipping (MAF).

- **Models** for tensor parameter distribution in amorphous materials.

  - Czjzek
  - Extended Czjzek

For more information, refer to the
[documentation](https://mrsimulator.readthedocs.io/en/stable/).

## Installation

    $ pip install mrsimulator

Please read our [installation document](https://mrsimulator.readthedocs.io/en/stable/installation.html) for details.

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/DeepanshS/mrsimulator-examples/master/test_file_v0.3.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://mrsimulator.readthedocs.io/en/master/_images/test_file.png)

## Reporting Bugs

The preferred location for submitting feature requests and bug reports is the [Github issue tracker](https://github.com/DeepanshS/mrsimulator/issues). Reports are also welcomed by directly contacting [Deepansh Srivastava](mailto:srivastava.89@osu.edu).

Discussions are welcome on [Github discussion](https://github.com/DeepanshS/mrsimulator/discussions)
