# -*- coding: utf-8 -*-
"""
Czjzek distribution (Shielding and Quadrupolar)
===============================================

In this example, we illustrate the simulation of spectrum originating from a Czjzek
distribution of traceless symmetric tensors. We show two cases, a Czjzek distribution
of the shielding and quadrupolar tensor parameters, respectively.

"""
# %%
# Import the required modules.
import matplotlib as mpl
import matplotlib.pyplot as plt
from mrsimulator import Simulator
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.models import czjzek_distribution
from mrsimulator.utils.collection import single_site_system_generator

# pre config the figures
mpl.rcParams["figure.figsize"] = [4.25, 3.0]

# %%
# Symmetric shielding tensor
# --------------------------
#
# Create the Czjzek distribution
# ''''''''''''''''''''''''''''''
#
# First, create a distribution of the zeta and eta parameters of the shielding tensors
# using the Czjzek distribution model as follows,
zeta_dist, eta_dist = czjzek_distribution(sigma=3.1415, n=5000)

# %%
# where `sigma` is the standard deviation of the second-rank tensor parameters, and `n`
# is the number of tensor parameters sampled from the Czjzek distribution model.
#
# The following is the plot of the Czjzek distribution.
plt.scatter(zeta_dist, eta_dist, s=2)
plt.xlabel(r"$\zeta$ / ppm")
plt.ylabel(r"$\eta$")
plt.ylim(0, 1)
plt.tight_layout()
plt.show()

# %%
# Simulate the spectrum
# '''''''''''''''''''''
#
# To quickly create the spin-systems from the above :math:`\zeta` and :math:`\eta`
# parameters, use the
# :func:`~mrsimulator.utils.collection.single_site_system_generator` utility function.
systems = single_site_system_generator(
    isotopes="13C", shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist}
)

# %%
# Here, the variable ``systems`` holds 5000 single-site spin systems.
print(len(systems))

# %%
# Create a simulator object and add the above system.
sim = Simulator()
sim.spin_systems = systems  # add the systems
sim.methods = [BlochDecaySpectrum(channels=["13C"])]  # add the method
sim.run()

# %%
# The following is the static spectrum arising from a Czjzek distribution of the
# second-rank traceless shielding tensors.
plt.figure(figsize=(4.25, 3.0))
ax = plt.gca(projection="csdm")
ax.plot(sim.methods[0].simulation, color="black")
plt.tight_layout()
plt.show()

# %%
# Quadrupolar tensor
# ------------------
#
# Create the Czjzek distribution
# ''''''''''''''''''''''''''''''
#
# Similarly, you may also create a Czjzek distribution of the electric field gradient
# (EFG) tensor parameters.
Cq_dist, eta_dist = czjzek_distribution(sigma=2.3, n=5000)

# %%
# Simulate the spectrum
# '''''''''''''''''''''
#
# Create the spin systems.
systems = single_site_system_generator(
    isotopes="71Ga", quadrupolar={"Cq": Cq_dist * 1e6, "eta": eta_dist}
)

# %%
# Create a simulator object and add the above system.
sim = Simulator()
sim.spin_systems = systems  # add the systems
sim.methods = [
    BlochDecayCentralTransitionSpectrum(
        channels=["71Ga"],
        magnetic_flux_density=4.8,  # in T
        spectral_dimensions=[{"count": 2048, "spectral_width": 1.2e6}],
    )
]  # add the method
sim.run()

# %%
# The following is the static spectrum arising from a Czjzek distribution of the
# second-rank traceless EFG tensors.
plt.figure(figsize=(4.25, 3.0))
ax = plt.gca(projection="csdm")
ax.plot(sim.methods[0].simulation, color="black")
ax.invert_xaxis()
plt.tight_layout()
plt.show()
