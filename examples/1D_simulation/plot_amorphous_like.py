#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amorphous-like materials
^^^^^^^^^^^^^^^^^^^^^^^^

29Si (I=1/2) simulation of amorphous-like material.
"""
# sphinx_gallery_thumbnail_number = 2
#%%
# In this section, we illustrate how ``Mrsimulator`` may be used in simulating
# amorphous materials. We do this by assuming a distribution of tensors. For
# example, consider the following,
import numpy as np
from scipy.stats import multivariate_normal

n = 4000
mean = [-100, 50, 0.15]  # given as [isotropic chemical shift in ppm, zeta in ppm, eta].
covariance = [[2.25, 0, 0], [0, 26.2, 0], [0, 0, 0.001]]  # same order as the mean.
iso, zeta, eta = multivariate_normal.rvs(mean=mean, cov=covariance, size=n).T

#%%
# Here, the coordinates ``iso``, ``zeta``, and ``eta`` are drawn
# from a three-dimension multivariate normal distribution of isotropic chemical
# shift, nuclear shielding anisotropy, and nuclear shielding asymmetry
# parameters, respectively. The mean of the distribution is given by the variable
# ``mean`` and holds a value of -100 ppm, 50 ppm, and 0.15 for the isotropic chemical
# shift, nuclear shielding anisotropy, and nuclear shielding asymmetry parameter,
# respectively. Similarly, the variable ``covariance`` holds the covariance matrix
# of the multivariate normal distribution. The two-dimensional plots from this
# three-dimensional distribution are shown below.

#%%
import matplotlib.pyplot as plt

_, ax = plt.subplots(1, 3, figsize=(9, 3))
ax[0].scatter(iso, zeta, color="black", s=0.5, alpha=0.3)
ax[0].set_xlabel("isotropic chemical shift / ppm")
ax[0].set_ylabel("shielding anisotropy / ppm")
ax[0].set_xlim(-120, -80)
ax[0].set_ylim(0, 100)

ax[1].scatter(iso, eta, color="black", s=0.5, alpha=0.3)
ax[1].set_xlabel("isotropic chemical shift / ppm")
ax[1].set_ylabel("shielding asymmetry")
ax[1].set_xlim(-120, -80)
ax[1].set_ylim(0, 1)

ax[2].scatter(zeta, eta, color="black", s=0.5, alpha=0.3)
ax[2].set_xlabel("shielding anisotropy / ppm")
ax[2].set_ylabel("shielding asymmetry")
ax[2].set_xlim(0, 100)
ax[2].set_ylim(0, 1)

plt.tight_layout()
plt.show()
#%%
#
# Let's create the site and isotopomer objects from these parameters.

#%%
from mrsimulator import Simulator, Site, Isotopomer, Dimension
from mrsimulator.methods import one_d_spectrum

isotopomers = []
for i, z, e in zip(iso, zeta, eta):
    site = Site(
        isotope="29Si",
        isotropic_chemical_shift=i,
        shielding_symmetric={"zeta": z, "eta": e},
    )
    isotopomers.append(Isotopomer(sites=[site]))

#%%
# Now, that we have the isotopomers, create a Simulator object and add the isotopomers.

#%%
sim = Simulator()
# add isotopomers
sim.isotopomers += isotopomers
# create and add a dimension
sim.dimensions += [
    Dimension(isotope="29Si", spectral_width=25000, reference_offset=-7000)
]

#%%
# Static line-shape
# -----------------
# Observe the static :math:`^{29}\text{Si}` line-shape simulation.

#%%
x, y = sim.run(method=one_d_spectrum)

#%%
# The plot of the corresponding spectrum.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# .. note::
#     The broad lineshape seen in the above spectrum is the result of
#     lineshapes arising from a distribution of tensors. In this case,
#     the lineshape is an integral of ``n`` individual spectra. There is no
#     lineshape broadening filter applied to the spectrum.

#%%
# Spinning sideband simulation at :math:`90^\circ`
# ------------------------------------------------
# Here is another example of a sideband simulation, spinning at a 90-degree angle,

#%%
sim.dimensions[0].rotor_frequency = 5000  # in Hz
sim.dimensions[0].rotor_angle = 1.57079  # 90 degree in radian
sim.config.number_of_sidebands = 8  # eight sidebands are sufficient for this example
x, y = sim.run(method=one_d_spectrum)

#%%
# and the corresponding plot.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# Spinning sideband simulation at magic angle
# -------------------------------------------

#%%
sim.dimensions[0].rotor_frequency = 1000  # in Hz
sim.dimensions[0].rotor_angle = 54.735 * np.pi / 180.0  # magic angle in radian
sim.config.number_of_sidebands = 16
x, y = sim.run(method=one_d_spectrum)

#%% and the corresponding plot.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()
