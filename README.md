# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://img.shields.io/pypi/v/mrsimulator.svg?style=flat&logo=pypi&logoColor=white)](https://pypi.python.org/pypi/mrsimulator) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mrsimulator)                                                                                                                                                                                                                                                                                                                                                                                                          |
| Build Status | [![GitHub Workflow Status](<https://img.shields.io/github/workflow/status/deepanshs/mrsimulator/CI?logo=GitHub>)](https://github.com/deepanshs/mrsimulator/actions) [![Read the Docs](https://img.shields.io/readthedocs/mrsimulator)](https://mrsimulator.readthedocs.io/en/stable/)                                                                                                                                                                                                                                                                                                                                         |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Metrics      | [![Total alerts](https://img.shields.io/lgtm/alerts/g/deepanshs/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/deepanshs/mrsimulator/alerts/) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/deepanshs/mrsimulator.svg?logo=lgtm)](https://lgtm.com/projects/g/deepanshs/mrsimulator/context:python) [![codecov](https://codecov.io/gh/deepanshs/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/deepanshs/mrsimulator) [![CodeFactor](https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge)](https://www.codefactor.io/repository/github/deepanshs/mrsimulator) |

`mrsimulator` is an open-source python package for fast computation/analysis of solid-state
magnetic resonance (NMR) spectra of both crystalline and amorphous materials.

**Why use mrsimulator?**

- It is open-source and free.
- It is a fast and versatile solid-state NMR simulator of one and two-dimensional static, MAS,
  and VAS spectra of nuclei experiencing chemical shift (nuclear shielding) and quadrupolar
  coupling interactions.
- It includes simulations of weakly coupled nuclei experiencing J and dipolar couplings.
- It is fully documented with a stable and simple API and is easily incorporated into your
  python scripts and web apps.
- It is compatible with modern python packages, such as scikit-learn, Keras, etc.
- Packages using mrsimulator:
  - [mrinversion](https://mrinversion.readthedocs.io/en/stable/)

## Install

```sh
pip install mrsimulator
```

Please refer to our [installation document](https://mrsimulator.readthedocs.io/en/latest/installation/users.html) for details.

#### A 1D static and MAS example

```py
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import BlochDecaySpectrum
import matplotlib.pyplot as plt

# Make Site and SpinSystem objects
H_site = Site(isotope="1H", shielding_symmetric={"zeta": 13.89, "eta": 0.25})
spin_system = SpinSystem(sites=[H_site])

# Make static and MAS one-pulse acquire Method objects
static = BlochDecaySpectrum(channels=["1H"]   )
mas = BlochDecaySpectrum(channels=["1H"], rotor_frequency=1000)  # in Hz

# Setup and run the Simulation object
sim = Simulator(spin_systems=[spin_system], methods=[static, mas])
sim.run()

# Plot the spectra
fig, ax = plt.subplots(1, 2, figsize=(6, 3), subplot_kw={"projection": "csdm"})
ax[0].plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax[0].set_title("Static")
ax[1].plot(sim.methods[1].simulation.real, color="black", linewidth=1)
ax[1].set_title("MAS")
plt.tight_layout()
plt.show()
```

This should produce the following figure.

<img src="docs/_static/test_file.png" width="800" />

<!-- ![alt text](docs/_static/test_file.png) -->

---

**Check out our extensive [documentation](https://mrsimulator.readthedocs.io/en/stable/index.html) and more example.**

[![Simulation](https://img.shields.io/badge/View-Simulation%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/examples/index.html)
[![Fitting](https://img.shields.io/badge/View-Fitting%20Examples-Purple?s=small)](https://mrsimulator.readthedocs.io/en/stable/fitting/index.html)

---

## Features

The `mrsimulator` package currently offers the following

- **Fast simulation** of one and two-dimensional solid-state NMR spectra.

- Simulation of **coupled and uncoupled spin system**

  - for spin I=1/2, and quadrupole I>1/2 nuclei
  - at arbitrary macroscopic magnetic flux density
  - at arbitrary rotor angles
  - at arbitrary spinning frequency

- A library of **NMR methods**,

  - 1D Bloch decay spectrum
  - 1D Bloch decay central transition spectrum
  - 2D Multi-Quantum Variable Angle Spinning (MQ-VAS)
  - 2D Satellite-transition Variable Angle Spinning (MQ-VAS)
  - 2D Dynamic Angle Spinning (DAS)
  - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT)
  - 2D Magic Angle Flipping (MAF)
  - Custom user-defined 1D and 2D methods (Method)

- **Models** for tensor parameter distribution in amorphous materials.

  - Czjzek
  - Extended Czjzek
  - Custom user-defined models

For more information, refer to the
[documentation](https://mrsimulator.readthedocs.io/en/stable/).

## Reporting Bugs

Submit bug reports or feature requests on the [Github issue tracker](https://github.com/deepanshs/mrsimulator/issues).

Discussions are welcome on the [Github discussion](https://github.com/deepanshs/mrsimulator/discussions) page.

## How to cite

If you use mrsimulator in your publication, please consider citing the following.

- Deepansh J. Srivastava, Matthew Giammar, Maxwell C. Venetos, Shyam Dwaraknath, Philip J. Grandinetti, & Alexis McCarthy. (2021). mrsimulator: v0.6.1. Zenodo. https://doi.org/10.5281/zenodo.5559730

- Srivastava DJ, Vosegaard T, Massiot D, Grandinetti PJ (2020) Core Scientific Dataset Model: A lightweight and portable model and file format for multi-dimensional scientific data. PLOS ONE 15(1): e0225953. https://doi.org/10.1371/journal.pone.0225953

_Additionally, if you use lmfit for least-squares fitting, consider citing the lmfit package._ Zenodo. http://doi.org/10.5281/zenodo.4516651
