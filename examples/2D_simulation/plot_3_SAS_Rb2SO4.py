#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rb2SO4, 87Rb (I=3/2) SAS
^^^^^^^^^^^^^^^^^^^^^^^^

87Rb (I=3/2) Switched-angle spinning (SAS) simulation.
"""
# %%
# The following is a switched-angle spinning (SAS) simulation of
# :math:`\text{Rb}_2\text{SO}_4` acquired at 9.4 T. There are two rubidium sites in
# :math:`\text{Rb}_2\text{SO}_4`. The tensor parameters for the two sites, used in this
# simulation, is taken from Shore `et. al.` [#f1]_.
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import MQVAS

# global plot configuration
font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 1

# %%
# **Step 1:** Create the site and the spin systems.
sites = [
    Site(
        isotope="87Rb",
        isotropic_chemical_shift=16,  # in ppm
        quadrupolar={"Cq": 5.3e6, "eta": 0.1},  # Cq in Hz
    ),
    Site(
        isotope="87Rb",
        isotropic_chemical_shift=40,  # in ppm
        quadrupolar={"Cq": 2.6e6, "eta": 1.0},  # Cq in Hz
    ),
]
spin_systems = [SpinSystem(sites=[s]) for s in sites]

# %%
# **Step 2:** Create a Multi Quantum variable angle spinning method. For SAS
# simulation, we use the same method as we used for MQ-MAS, but with modified
# parameters, as shown below.
method = MQVAS(
    channels=["87Rb"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=[
        {
            "count": 256,
            "spectral_width": 3.2e4,  # in Hz
            "reference_offset": -1e3,  # in Hz
            "label": "90 dimension",
            "events": [{"rotor_angle": 90 * np.pi / 180}],  # in radians
        },
        {
            "count": 256,
            "spectral_width": 22e3,  # in Hz
            "reference_offset": -4e3,  # in Hz
            "label": "MAS dimension",
            "events": [{"rotor_angle": 54.74 * np.pi / 180}],  # in radians
        },
    ],
)

# %%
# **Step 3:** Create the Simulator object, add the method and spin system objects, and
# run the simulation
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method.

# configure the simulator object
sim.run()

# %%
# **Step 4:**
# The plot of the simulation.
data = sim.methods[0].simulation
ax = plt.subplot(projection="csdm")
ax.imshow(data / data.max(), aspect="auto", cmap="gist_ncar_r")
ax.invert_xaxis()
ax.set_ylim(-70, 90)
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Shore, J.S., Wang, S.H., Taylor, R.E., Bell, A.T., Pines, A. Determination of
#       quadrupolar and chemical shielding tensors using solid-state two-dimensional NMR
#       spectroscopy, J. Chem. Phys. (1996)  **105** *21*, 9412.
#       `DOI: 10.1063/1.472776 <https://doi.org/10.1063/1.472776>`_
