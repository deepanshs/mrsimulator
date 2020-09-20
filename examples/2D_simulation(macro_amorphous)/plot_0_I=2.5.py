#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amorphous material, 27Al (I=5/2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

27Al (I=5/2) simulation of amorphous material.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mrsimulator import Simulator
from mrsimulator.methods import ThreeQ_VAS
from mrsimulator.utils.collection import single_site_system_generator
from scipy.stats import multivariate_normal

# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]

# %%
# In this section, we illustrate the simulation of a quadrupolar spectrum arising from
# a distribution of the electric field gradient (EFG) tensors from an amorphous
# material. We proceed by assuming a multi-variate normal distribution, as follows,

mean = [58, 4, 0.98]  # given as [isotropic chemical shift in ppm, Cq in MHz, eta].
covariance = [[6, 0, 0], [0, 1.4, 0], [0, 0, 0.35]]  # same order as the mean.

# range of coordinates along the three dimensions
iso_range = np.arange(100) / 1.5 + 30  # in ppm
Cq_range = np.arange(100) / 4 - 5  # in MHz
eta_range = np.arange(7) / 6

# The coordinates grid
iso, Cq, eta = np.meshgrid(iso_range, Cq_range, eta_range, indexing="ij")
pos = np.asarray([iso, Cq, eta]).T

# Three-dimensional probability distribution function.
pdf = multivariate_normal(mean=mean, cov=covariance).pdf(pos).T

# %%
# Here, ``iso``, ``Cq``, and ``eta`` are the isotropic chemical shift, the quadrupolar
# coupling constant, and quadrupolar asymmetry coordinates of the 3D-grid
# system over which the multivariate normal probability distribution is evaluated. The
# mean of the distribution is given by the variable ``mean`` and holds a value of 20
# ppm, 6.5 MHz, and 0.3 for the isotropic chemical shift, the quadrupolar coupling
# constant, and quadrupolar asymmetry parameter, respectively. Similarly, the variable
# ``covariance`` holds the covariance matrix of the multivariate normal distribution.
# The two-dimensional projections from this three-dimensional distribution are shown
# below.
_, ax = plt.subplots(1, 3, figsize=(9, 3))

# isotropic shift v.s. quadrupolar coupling constant
ax[0].contourf(Cq_range, iso_range, pdf.sum(axis=2))
ax[0].set_xlabel("Cq / MHz")
ax[0].set_ylabel("isotropic chemical shift / ppm")

# isotropic shift v.s. quadrupolar asymmetry
ax[1].contourf(eta_range, iso_range, pdf.sum(axis=1))
ax[1].set_xlabel(r"quadrupolar asymmetry, $\eta$")
ax[1].set_ylabel("isotropic chemical shift / ppm")

# quadrupolar coupling constant v.s. quadrupolar asymmetry
ax[2].contourf(eta_range, Cq_range, pdf.sum(axis=0))
ax[2].set_xlabel(r"quadrupolar asymmetry, $\eta$")
ax[2].set_ylabel("Cq / MHz")

plt.tight_layout()
plt.show()

# %%
# Let's create the site and spin system objects from these parameters. Note, we create
# single-site spin systems for optimum performance.
# Use the :func:`~mrsimulator.utils.collection.single_site_system_generator` utility
# function to generate single-site spin systems.
spin_systems = single_site_system_generator(
    isotopes="27Al",
    isotropic_chemical_shifts=iso,
    quadrupolar={"Cq": Cq * 1e6, "eta": eta},  # Cq in Hz
    abundance=pdf,
)
len(spin_systems)

# %%
# Static spectrum
# ---------------
# Observe the static :math:`^{27}\text{Al}` NMR spectrum simulation. First,
# create a central transition selective Bloch decay spectrum method.
mqvas = ThreeQ_VAS(
    channels=["27Al"],
    spectral_dimensions=[
        {
            "count": 512,
            "spectral_width": 26718.475776,  # in Hz
            "reference_offset": -4174.76184,  # in Hz
            "label": "Isotropic dimension",
        },
        {
            "count": 512,
            "spectral_width": 2e4,  # in Hz
            "reference_offset": 2e3,  # in Hz
            "label": "MAS dimension",
        },
    ],
)

# %%
# Create the simulator object and add the spin systems and method.
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [mqvas]  # add the method
sim.config.number_of_sidebands = 1
sim.run()

data = sim.methods[0].simulation

# %%
# The plot of the corresponding spectrum.
ax = plt.subplot(projection="csdm")
ax.imshow(data / data.max(), cmap="gist_ncar_r", aspect="auto")
ax.set_ylim(-20, -50)
ax.set_xlim(80, 20)
plt.tight_layout()
plt.show()
