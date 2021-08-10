# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://img.shields.io/pypi/v/mrsimulator.svg?style=flat&logo=pypi&logoColor=white)](https://pypi.python.org/pypi/mrsimulator) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mrsimulator)                                                                                                                                                                                                                                                                                                                                                                                                          |
| Build Status | [![GitHub Workflow Status](<https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI?logo=GitHub>)](https://github.com/deepanshs/mrsimulator/actions) [![Read the Docs](https://img.shields.io/readthedocs/mrsimulator)](https://mrsimulator.readthedocs.io/en/stable/)                                                                                                                                                                                                                                                                                                                                         |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Metrics      | [![Total alerts](https://img.shields.io/lgtm/alerts/g/deepanshs/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/deepanshs/mrsimulator/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/deepanshs/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/deepanshs/mrsimulator/context:python) [![codecov](https://codecov.io/gh/deepanshs/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/deepanshs/mrsimulator) [![CodeFactor](https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge)](https://www.codefactor.io/repository/github/deepanshs/mrsimulator) |

`mrsimulator` is an open-source python package for fast computation/analysis of solid-state
magnetic resonance (NMR) spectra of both crystalline and amorphous materials. The core
of the `mrsimulator` library is written in C, wrapped, and made available in python.

---

**Check out our extensive [documentation](https://mrsimulator.readthedocs.io/en/stable/index.html).**

[![Simulation](https://img.shields.io/badge/View-Simulation%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/examples/index.html)
[![Fitting](https://img.shields.io/badge/View-Fitting%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/fitting/index.html)

---

:warning: The package is currently under development. We advise using with caution. Bug report are greatly appreciated.

## Why use mrsimulator?

- It is open-source and free.
- It is a fast and versatile solid-state NMR simulator of one-dimensional static, MAS,
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- It include simulations of weakly coupled nuclei experiencing J and dipolar couplings.
- It is fully documented with a stable and simple API and is easily incorporated into your
  python scripts and web apps.
- It is compatible with modern python packages, such as scikit-learn, Keras, etc.
- Packages using mrsimulator -
  - [mrinversion](https://mrinversion.readthedocs.io/en/stable/)

### Features

The `mrsimulator` package currently offers the following

- **Fast simulation** of one-dimensional solid-state NMR spectra. See our
  [benchmark results](https://mrsimulator.readthedocs.io/en/stable/benchmark.html).

- Simulation of **coupled and uncoupled spin system**

  - for spin I=1/2, and quadrupole I>1/2 nuclei,
  - at arbitrary macroscopic magnetic flux density,
  - at arbitrary rotor angles, and
  - at arbitrary spinning frequency.

- A library of **NMR methods**,

  - 1D Bloch decay spectrum,
  - 1D Bloch decay central transition spectrum,
  - 2D Multi-Quantum Variable Angle Spinning (MQ-VAS),
  - 2D Satellite-transition Variable Angle Spinning (MQ-VAS),
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

Please read our [installation document](https://mrsimulator.readthedocs.io/en/latest/installation/users.html) for details.

## Check your build

If the installation is successful, you should be able to run the following
[test file](https://raw.github.com/deepanshs/mrsimulator-examples/master/test_file_v0.3.py?raw=true)
in your terminal.

    $ python test_file.py

This should produce the following figure.

![alt text](https://mrsimulator.readthedocs.io/en/master/_images/test_file.png)

## Reporting Bugs

The preferred location for submitting feature requests and bug reports is the [Github issue tracker](https://github.com/deepanshs/mrsimulator/issues). Reports are also welcomed by directly contacting [Deepansh Srivastava](mailto:srivastava.89@osu.edu).

Discussions are welcome on [Github discussion](https://github.com/deepanshs/mrsimulator/discussions)

## How to cite

If you use mrsimulator in your publication, please consider citing the following.

- Deepansh J. Srivastava, Maxwell Venetos, Philip J. Grandinetti, Shyam Dwaraknath, & Alexis McCarthy. (2021, May 26). mrsimulator: v0.6.0 (Version v0.6.0). Zenodo. http://doi.org/10.5281/zenodo.4814638

- Srivastava DJ, Vosegaard T, Massiot D, Grandinetti PJ (2020) Core Scientific Dataset Model: A lightweight and portable model and file format for multi-dimensional scientific data. PLOS ONE 15(1): e0225953. https://doi.org/10.1371/journal.pone.0225953

_Additionally, if you use lmfit for least-squares fitting, consider citing the lmfit package._

- Matt Newville; Renee Otten; Andrew Nelson; Antonino Ingargiola; Till Stensitzki; Dan Allan; Austin Fox; Faustin Carter; Michał; Dima Pustakhod; lneuhaus; Sebastian Weigand; Ray Osborn; Glenn; Christoph Deil; Mark; Allan L. R. Hansen; Gustavo Pasquevich; Leon Foks; Nicholas Zobrist; Oliver Frost; Alexandre Beelen; Stuermer; kwertyops; Anthony Polloreno; Shane Caldwell; Anthony Almarza; Arun Persaud; Ben Gamari; Benjamin F. Maier. (2021, February 7). lmfit/lmfit-py 1.0.2 (Version 1.0.2). Zenodo. http://doi.org/10.5281/zenodo.4516651
