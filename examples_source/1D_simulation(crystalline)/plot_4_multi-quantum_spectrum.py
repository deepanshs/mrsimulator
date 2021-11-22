#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Arbitrary spin transition (multi-quantum)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

³³S (I=5/2) quadrupolar spectrum simulation.
"""
# %%
# Simulate a triple quantum spectrum.
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import Method1D
from mrsimulator.method.event import SpectralEvent
from mrsimulator.spin_system.tensors import SymmetricTensor

# sphinx_gallery_thumbnail_number = 2

# %%
# Create a single-site arbitrary spin system.
site = Site(
    name="27Al",
    isotope="27Al",
    isotropic_chemical_shift=35.7,  # in ppm
    quadrupolar=SymmetricTensor(Cq=2.959e6, eta=0.98),  # Cq is in Hz
)
spin_system = SpinSystem(sites=[site])

# %%
# Selecting the triple-quantum transition
# ---------------------------------------
#
# For single-site spin-5/2 spin system, there are three triple-quantum transition
#
# - :math:`|1/2\rangle\rightarrow|-5/2\rangle` (:math:`P=-3, D=6`)
# - :math:`|3/2\rangle\rightarrow|-3/2\rangle` (:math:`P=-3, D=0`)
# - :math:`|5/2\rangle\rightarrow|-1/2\rangle` (:math:`P=-3, D=-6`)
#
# To select one or more triple-quantum transitions, assign the respective value of P and
# D to the `transition_query`. Here, we select the symmetric triple-quantum transition.
method = Method1D(
    name="Arbitrary Transition Method",
    channels=["27Al"],
    magnetic_flux_density=21.14,  # in T
    rotor_frequency=1e9,  # in Hz
    spectral_dimensions=[
        dict(
            count=1024,
            spectral_width=5e3,  # in Hz
            reference_offset=2.5e4,  # in Hz
            events=[
                SpectralEvent(
                    # symmetric triple quantum transitions
                    transition_query=[dict(P=[-3], D=[0])]
                ),
            ],
        )
    ],
)

# A graphical representation of the method object.
plt.figure(figsize=(5, 3))
method.plot()
plt.show()

# %%
# Create the Simulator object and add the method and the spin system object.
sim = Simulator()
sim.spin_systems = [spin_system]  # add the spin system
sim.methods = [method]  # add the method
sim.run()

# The plot of the simulation before signal processing.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
