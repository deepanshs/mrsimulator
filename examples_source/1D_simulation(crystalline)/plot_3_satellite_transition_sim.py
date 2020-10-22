#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simulate arbitrary transitions (single-quantum)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

27Al (I=5/2) quadrupolar spectrum simulation.
"""
# %%
# The ``mrsimulator`` library does not offer any pre-defined method for simulating
# individual transitions. A BlochDecaySpectrum method simulates all single quantum
# transitions, while a BlochDecayCentralTransitionSpectrum method only simulates the
# central transition. In this example, we show how you can simulate
# any arbitrary transition using the generic Method1D method.
import matplotlib as mpl
import matplotlib.pyplot as plt
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import Method1D

# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
# sphinx_gallery_thumbnail_number = 2

# %%
# Create a single-site arbitrary spin system.
site = Site(
    name="27Al",
    isotope="27Al",
    isotropic_chemical_shift=35.7,  # in ppm
    quadrupolar={"Cq": 5.959e6, "eta": 0.32},  # Cq is in Hz
)
spin_system = SpinSystem(sites=[site])

# %%
# Selecting the inner-satellite transition
# ----------------------------------------
#
# The arguments of the following Method1D object is the same as the BlochDecaySpectrum
# method. One extra argument is the `events` item in the `spectral_dimension` object.
# The event object is where you define the `transition_query` to select one or more
# transitions to simulate. The two attributes of the `transition_query` are `P` and `D`,
# which are given as, :math:`m_f-m_i` and :math:`m_f^2 - m_i^2`, where :math:`m_f` and
# :math:`m_i` are the spin quantum numbers for the final and initial energy states.
#
# In the following example, we assign the values of P and D as -1 and 2, respectively.
# In the case of a single-site spin 5/2 spin system, there is only one transition,
# :math:`|-1/2\rangle\rightarrow|-3/2\rangle`, that satisfy this query selection
# criterion and thus will be selected.
#
# For illustrative purposes, let's look at the infinite speed spectrum from this
# satellite transition.
method = Method1D(
    channels=["27Al"],
    magnetic_flux_density=21.14,  # in T
    rotor_frequency=1e9,  # in Hz
    spectral_dimensions=[
        {
            "count": 1024,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": 1e4,  # in Hz
            "events": [
                {"transition_query": {"P": [-1], "D": [2]}}  # <-- select transitions
            ],
        }
    ],
)

# %%
# Create the Simulator object and add the method and the spin system object.
sim = Simulator()
sim.spin_systems += [spin_system]  # add the spin system
sim.methods += [method]  # add the method

# %%
# Simulate the spectrum.
sim.run()

# The plot of the simulation before signal processing.
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Selecting both inner and outer-satellite transition
# ---------------------------------------------------
# You may use the same transition query selection criterion to select multiple
# transitions. Consider the following transitions with respective P and D values.
#
# - :math:`|-1/2\rangle\rightarrow|-3/2\rangle` (:math:`P=-1, D=2`)
# - :math:`|-3/2\rangle\rightarrow|-5/2\rangle` (:math:`P=-1, D=4`)
method2 = Method1D(
    channels=["27Al"],
    magnetic_flux_density=21.14,  # in T
    rotor_frequency=1e9,  # in Hz
    spectral_dimensions=[
        {
            "count": 1024,
            "spectral_width": 1e4,  # in Hz
            "reference_offset": 1e4,  # in Hz
            "events": [
                {"transition_query": {"P": [-1], "D": [2, 4]}}  # <-- select transitions
            ],
        }
    ],
)

# %%
# Update the method object in the Simulator object.
sim.methods[0] = method2  # add the method

# %%
# Simulate the spectrum.
sim.run()

# The plot of the simulation before signal processing.
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
