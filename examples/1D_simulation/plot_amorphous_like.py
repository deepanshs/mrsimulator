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

n = 1000
iso = np.random.normal(loc=-100.0, scale=5.0, size=n)
zeta = np.random.normal(loc=50.0, scale=12.12, size=n)
eta = np.random.normal(loc=0.5, scale=0.1, size=n)

#%%
# Here, we have created three Gaussian distributions, one for each isotropic
# chemical shift, ``iso``, shielding anisotropy, ``zeta``, and shielding
# asymmetry, ``eta``. The isotropic chemical shift distribution is centered at
# -100 ppm with a standard deviation of 5 ppm. Similarly, the ``zeta`` and
# ``eta`` distributions are centered at 50 ppm and 0.5 with 12.12 ppm and 0.1
# standard deviations, respectively. A total of 1000 points is sampled from the
# three distributions, resulting in an uncorrelated ``iso``-``zeta``-``eta``
# distribution, shown below

#%%
import matplotlib.pyplot as plt

_, ax = plt.subplots(1, 3, figsize=(9, 3))
ax[0].scatter(iso, zeta, color="black", s=0.5)
ax[0].set_xlabel("isotropic chemical shift / ppm")
ax[0].set_ylabel("shielding anisotropy / ppm")
ax[0].set_xlim(-120, -80)
ax[0].set_ylim(0, 100)

ax[1].scatter(iso, eta, color="black", s=0.5)
ax[1].set_xlabel("isotropic chemical shift / ppm")
ax[1].set_ylabel("shielding asymmetry")
ax[1].set_xlim(-120, -80)
ax[1].set_ylim(0, 1)

ax[2].scatter(zeta, eta, color="black", s=0.5)
ax[2].set_xlabel("shielding anisotropy / ppm")
ax[2].set_ylabel("shielding asymmetry")
ax[2].set_xlim(0, 100)
ax[2].set_ylim(0, 1)

plt.tight_layout()
plt.show()
#%%
#
# Let's create the site and isotopomers objects from these parameters.

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
# Now, that we have a 1000 isotopomers, let's create the Simulator object and add
# these isotopomers.

#%%
sim = Simulator()
# add isotopomers
sim.isotopomers = isotopomers
# create and add a dimension
sim.dimensions = [
    Dimension(isotope="29Si", spectral_width=25000, reference_offset=-7000)
]

#%%
# Static line-shape
# -----------------
# Observe the static line-shape simulation.

#%%
x, y = sim.run(method=one_d_spectrum)

#%%
# The plot of the corresponding spectrum.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.tight_layout()
plt.show()

#%%
# .. note::
#     The broad lineshape seen in the above spectrum is the result of
#     lineshapes arising from a distribution of tensors. In this case,
#     the lineshape is an integral of 1000 individual spectrum. There is no
#     lineshape broadening filter applied to the spectrum.

#%%
# Spinning sidebands
# ------------------
# Here is another example of a sideband simulation, spinning at a 90-degree angle,

#%%
sim.dimensions[0].rotor_frequency = 5000  # in Hz
sim.dimensions[0].rotor_angle = 1.57079  # 90 degree in radian
sim.config.number_of_sidebands = 8  # eight sidebands are sufficient for this example
x, y = sim.run(method=one_d_spectrum)

#%% and the corresponding plot.

#%%
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.tight_layout()
plt.show()
