#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amorphous-like materials
^^^^^^^^^^^^^^^^^^^^^^^^

29Si (I=1/2) simulation of amorphous-like material.
"""
# sphinx_gallery_thumbnail_number = 2
#%%
# One of the advantages of ``Mrsimulator`` is that it is a very fast NMR line-shape
# simulation library. We can exploit this feature to simulate bulk line-shapes and
# eventually model amorphous materials.
#
# In this section, we illustrate how ``Mrsimulator`` may be used in simulating
# NMR spectrum of amorphous materials. We accomplish this by assuming a distribution
# of interaction tensors in amorphous materials. For example, consider a tri-variate
# normal distribution of shielding tensor parameters, `i.e.`, the isotropic chemical
# shift, the anisotropy parameter, :math:`\zeta`, and the asymmetry parameter,
# :math:`\eta`, as follows,
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
# Let's create the sites and isotopomer objects from these parameters.

#%%
from mrsimulator import Simulator, Site, Isotopomer

isotopomers = []
for i, z, e in zip(iso, zeta, eta):
    site = Site(
        isotope="29Si",
        isotropic_chemical_shift=i,
        shielding_symmetric={"zeta": z, "eta": e},
    )
    isotopomers.append(Isotopomer(sites=[site]))

#%%
# Let's also create the Bloch decay spectrum method.
from mrsimulator.methods import BlochDecaySpectrum

method = BlochDecaySpectrum(
    channels=["29Si"],
    spectral_dimensions=[
        {"spectral_width": 25000, "reference_offset": -7000}
    ],  # values in Hz
)

#%%
# The above method simulates static :math:`^{29}\text{Si}` line-shapes at 9.4 T field
# (default value).
#
# Now, that we have the isotopomers and the method, create the simulator object and
# all the respective objects.

sim = Simulator()
sim.isotopomers += isotopomers  # add isotopomers
sim.methods += [method]  # add method

#%%
# Static line-shape
# -----------------
# Observe the static :math:`^{29}\text{Si}` line-shape simulation.

#%%
sim.run()

#%%
# The plot of the simulation.

x, y = sim.methods[0].simulation.to_list()
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
#     the lineshape is an integral of ``n=4000`` individual spectra. There is no
#     lineshape broadening filter applied to the spectrum.

#%%
# Spinning sideband simulation at :math:`90^\circ`
# ------------------------------------------------
# Here is another example of a sideband simulation, spinning at a 90-degree angle.

sim.methods[0] = BlochDecaySpectrum(
    channels=["29Si"],
    rotor_frequency=5000,  # in Hz
    rotor_angle=1.57079,  # in rads, equivalent to 90 deg
    spectral_dimensions=[
        {"spectral_width": 25000, "reference_offset": -7000}  # values in Hz
    ],
)
sim.config.number_of_sidebands = 8  # eight sidebands are sufficient for this example
sim.run()

#%%
# The plot of the simulation.

x, y = sim.methods[0].simulation.to_list()
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

sim.methods[0] = BlochDecaySpectrum(
    channels=["29Si"],
    rotor_frequency=1000,  # in Hz
    rotor_angle=54.735 * np.pi / 180.0,  # in rads
    spectral_dimensions=[
        {"spectral_width": 25000, "reference_offset": -7000}  # values in Hz
    ],
)
sim.config.number_of_sidebands = 16  # sixteen sidebands are sufficient for this example
sim.run()

#%%
# The plot of the simulation.

x, y = sim.methods[0].simulation.to_list()
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()
