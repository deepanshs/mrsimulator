#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Potassium Sulfate
^^^^^^^^^^^^^^^^^

33S (I=3/2) quadrupolar line-shape simulation.
"""
# %%
# The following example is the :math:`^{33}\text{S}` NMR line-shape simulation of
# potassium sulfate (:math:`\text{K}_2\text{SO}_4`). The quadrupole tensor parameters
# for :math:`^{33}\text{S}` is obtained from Moudrakovski et. al. [#f3]_
import matplotlib.pyplot as plt
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import one_d_spectrum

#%%
# **Step 1** Create sites, in this case, just the one.

S33 = Site(
    name="33S",
    isotope="33S",
    isotropic_chemical_shift=335.7,
    quadrupolar={"Cq": 0.959e6, "eta": 0.42},
)

#%%
# **Step 2** Create isotopomers from these sites.

isotopomer = Isotopomer(sites=[S33])

#%%
# **Step 3** Create the dimension.


dimension = Dimension(
    isotope="33S",
    magnetic_flux_density=21.14,  # in T
    number_of_points=2046,
    spectral_width=5000,  # in Hz
    reference_offset=22500,  # in Hz
    rotor_frequency=14000,  # in Hz
)

#%%
# **Step 4** Create the Simulator object and add dimension and isotopomer objects.

sim_K2SO3 = Simulator()

# add isotopomers
sim_K2SO3.isotopomers += [isotopomer]

# add dimensions
sim_K2SO3.dimensions += [dimension]

#%%
# **Step 5** Simulate the spectrum.

x, y = sim_K2SO3.run(method=one_d_spectrum)

#%%
# **Step 6** Plot.

plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# .. [#f3] Moudrakovski, I., Lang, S., Patchkovskii, S., and Ripmeester, J. High field
#       :math:`^{33}\text{S}` solid state NMR and first-principles calculations in
#       potassium sulfates. J. Phys. Chem. A, 2010, **114**, *1*, 309–316.
#       `DOI: 10.1021/jp908206c <https://doi.org/10.1021/jp908206c>`_
