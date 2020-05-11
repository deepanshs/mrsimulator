#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amorphous-like materials (quadrupolar)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

27Al (I=5/2) simulation of amorphous-like material.
"""
# sphinx_gallery_thumbnail_number = 2
#%%
# In this section, we illustrate the simulation of a quadrupolar spectrum
# arising from amorphous materials. We proceed by assuming a
# distribution of electric field gradient (EFG) tensors, as follows,
import numpy as np
from scipy.stats import multivariate_normal

n = 4000
mean = [20, 6.5, 0.3]  # given as [isotropic chemical shift in ppm, Cq in MHz, eta].
covariance = [[1.98, 0, 0], [0, 4.9, 0], [0, 0, 0.0016]]  # same order as the mean.
iso, Cq, eta = multivariate_normal.rvs(mean=mean, cov=covariance, size=n).T

#%%
# Here, the coordinates ``iso``, ``Cq``, and ``eta`` are drawn
# from a three-dimension multivariate normal distribution of isotropic chemical
# shift, electric quadrupole coupling constant, and quadrupole asymmetry
# parameters, respectively. The mean of the distribution is given by the variable
# ``mean`` and holds a value of 20 ppm, 6.5 MHz, and 0.3 for the isotropic chemical
# shift, electric quadrupole coupling constant, and quadrupole asymmetry parameter,
# respectively. Similarly, the variable ``covariance`` holds the covariance matrix
# of the multivariate normal distribution. The two-dimensional plots from this
# three-dimensional distribution are shown below.

#%%
import matplotlib.pyplot as plt

_, ax = plt.subplots(1, 3, figsize=(9, 3))
ax[0].scatter(iso, Cq, color="black", s=0.5, alpha=0.3)
ax[0].set_xlabel("isotropic chemical shift / ppm")
ax[0].set_ylabel("Cq / MHz")
ax[0].set_xlim(10, 30)
ax[0].set_ylim(-10, 20)

ax[1].scatter(iso, eta, color="black", s=0.5, alpha=0.3)
ax[1].set_xlabel("isotropic chemical shift / ppm")
ax[1].set_ylabel("quadrupolar asymmetry")
ax[1].set_xlim(10, 30)
ax[1].set_ylim(0, 1)

ax[2].scatter(Cq, eta, color="black", s=0.5, alpha=0.3)
ax[2].set_xlabel("Cq / MHz")
ax[2].set_ylabel("quadrupolar asymmetry")
ax[2].set_xlim(-10, 20)
ax[2].set_ylim(0, 1)

plt.tight_layout()
plt.show()
#%%
#
# Let's create the site and isotopomer objects from these parameters.

#%%
from mrsimulator import Simulator, Site, Isotopomer

isotopomers = []
for i, c, e in zip(iso, Cq, eta):
    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=i,
        quadrupolar={"Cq": c * 1e6, "eta": e},  # Cq in Hz
    )
    isotopomers.append(Isotopomer(sites=[site]))

#%%
# Now, that we have the isotopomers, create a Simulator object and add the isotopomers.

#%%
sim = Simulator()
# add isotopomers
sim.isotopomers += isotopomers
# create and add a central transition selective Bloch decay spectrum method.
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum

method = BlochDecayCentralTransitionSpectrum(
    channels=["27Al"], dimensions=[{"spectral_width": 80000}]
)

sim.methods += [method]

#%%
# Static line-shape
# -----------------
# Observe the static :math:`^{27}\text{Al}` line-shape simulation.

#%%
sim.run()
sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")
x, y = sim.methods[0].simulation.to_list()

#%%
# The plot of the corresponding spectrum.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("$^{27}$Al frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# Spinning sideband simulation at magic angle
# -------------------------------------------
# Simulation of the same isotopomer system at magic angle and spinning at 25 kHz.

#%%
sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency = 25000  # in Hz
sim.methods[0].spectral_dimensions[0].spectral_width = 30000  # in Hz
sim.methods[0].spectral_dimensions[0].reference_offset = -4000  # in Hz
sim.methods[0].spectral_dimensions[0].events[0].rotor_angle = (
    54.735 * np.pi / 180.0
)  # magic angle in radian
sim.config.number_of_sidebands = 4
sim.run()
sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")
x, y = sim.methods[0].simulation.to_list()

#%% and the corresponding plot.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("$^{27}$Al frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()
