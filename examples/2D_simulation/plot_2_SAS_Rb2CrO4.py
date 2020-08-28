#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rb2CrO4, 87Rb (I=3/2) SAS
^^^^^^^^^^^^^^^^^^^^^^^^^

87Rb (I=3/2) Switched-angle spinning (SAS) simulation.
"""
# %%
# The following is a switched-angle spinning (SAS) simulation of
# :math:`\text{Rb}_2\text{CrO}_4`. While :math:`\text{Rb}_2\text{CrO}_4` has two
# rubidium sites, the site with the smaller quadrupolar interaction was selectively
# observed and reported by Shore `et. al.` [#f1]_. The following is the simulation
# of the reported tensor parameters.
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
site = Site(
    isotope="87Rb",
    isotropic_chemical_shift=-7,  # in ppm
    shielding_symmetric={"zeta": 110, "eta": 0},
    quadrupolar={
        "Cq": 3.5e6,  # in Hz
        "eta": 0.3,
        "alpha": 0,  # in rads
        "beta": 70 * np.pi / 180,  # in rads
        "gamma": 0,  # in rads
    },
)
spin_system = SpinSystem(sites=[site])

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
            "spectral_width": 1.2e4,  # in Hz
            "reference_offset": -3e3,  # in Hz
            "label": "70.12 dimension - 1",
            "events": [
                {
                    "rotor_angle": 70.12 * np.pi / 180,  # in radians
                    "transition_query": {"P": [-1], "D": [0]},
                }
            ],
        },
        {
            "count": 512,
            "spectral_width": 8e3,  # in Hz
            "reference_offset": -4e3,  # in Hz
            "label": "MAS dimension - 0",
            "events": [
                {
                    "rotor_angle": 54.74 * np.pi / 180,  # in radians
                    "transition_query": {"P": [-1], "D": [0]},
                }
            ],
        },
    ],
)

# %%
# **Step 3:** Create the Simulator object, add the method and spin system objects, and
# run the simulation
sim = Simulator()
sim.spin_systems += [spin_system]  # add the spin systems
sim.methods += [method]  # add the method.

# configure the simulator object
sim.config.integration_volume = "hemisphere"
sim.config.number_of_sidebands = 1
sim.run()

# %%
# **Step 4:**
# The plot of the simulation.
data = sim.methods[0].simulation
ax = plt.subplot(projection="csdm")
ax.imshow(data / data.max(), aspect="auto", cmap="gist_ncar_r")
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Shore, J.S., Wang, S.H., Taylor, R.E., Bell, A.T., Pines, A. Determination of
#       quadrupolar and chemical shielding tensors using solid-state two-dimensional NMR
#       spectroscopy, J. Chem. Phys. (1996)  **105** *21*, 9412.
#       `DOI: 10.1063/1.472776 <https://doi.org/10.1063/1.472776>`_
