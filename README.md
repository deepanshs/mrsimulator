# The Mrsimulator project

|              |                                                                                                                                                                                                                                                                                                                                                                            |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Deployment   | [![PyPI version](https://img.shields.io/pypi/v/mrsimulator.svg?style=flat&logo=pypi&logoColor=white)](https://pypi.python.org/pypi/mrsimulator) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mrsimulator)                                                                                                                                               |
| Build Status | [![CI](https://github.com/deepanshs/mrsimulator/actions/workflows/continuous-integration-pip.yml/badge.svg?branch=master)](https://github.com/deepanshs/mrsimulator/actions/workflows/continuous-integration-pip.yml) [![Read the Docs](https://img.shields.io/readthedocs/mrsimulator)](https://mrsimulator.readthedocs.io/en/latest/) |
| License      | [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)                                                                                                                                                                                                                                                  |
| Metrics      | [![codecov](https://codecov.io/gh/deepanshs/mrsimulator/branch/master/graph/badge.svg)](https://codecov.io/gh/deepanshs/mrsimulator) [![CodeFactor](https://www.codefactor.io/repository/github/deepanshs/mrsimulator/badge)](https://www.codefactor.io/repository/github/deepanshs/mrsimulator)                                                                           |

Shortly after the birth of Nuclear Magnetic Resonance (NMR) spectroscopy, it was realized that spin
and spatial degrees of freedom could be manipulated on a time scale faster than the coherence
lifetimes of the nuclear spin transitions. This led to an explosion of multi-pulse and sample
reorientation methodologies in magnetic resonance for probing the structure and dynamics of matter
over a wide range of length and time scales.

Numerical simulations of the NMR spectra from these methods have long been a critical
part of their analyses. The most robust and rigorous numerical approaches employ the full density
operator, ideal for dealing with finite pulse effects, weak to intermediate to strong couplings,
non-commuting Hamiltonians, and relaxation and exchange processes. However, such approaches can be
highly inefficient, particularly when Hamiltonians commute, pulses are ideal, and transverse relaxation
can be treated as an ad-hoc line broadening. `mrsimulator`, an open-source python package, achieves
high benchmarks in spectral simulations and analyses by limiting itself to these simpler situations.
Fortunately, working within this limit only prevents `mrsimulator` from modeling spectra of a small
fraction of popular NMR methods. The efficiency gains with this approach over conventional density
operator simulations are tremendous.

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

Please refer to our [installation documentation](https://mrsimulator.readthedocs.io/en/stable/installation/users.html) for details.

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
  - 2D isotropic/anisotropic sideband correlation spectrum (e.g. PASS and MAT)
  - 2D Magic Angle Flipping (MAF)
  - 2D Dynamic Angle Spinning (DAS)
  - Custom user-defined methods (Method)

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

- Please use the GitHub citation tool to cite this repository. The tool in located in the About section under the `Cite this repository` category.

- Srivastava DJ, Vosegaard T, Massiot D, Grandinetti PJ (2020) Core Scientific Dataset Model: A lightweight and portable model and file format for multi-dimensional scientific dataset. PLOS ONE 15(1): e0225953. https://doi.org/10.1371/journal.pone.0225953

_Additionally, if you use lmfit for least-squares fitting, consider citing the lmfit package._ Zenodo. https://doi.org/10.5281/zenodo.4516651
