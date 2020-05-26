#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coesite
^^^^^^^

17O (I=5/2) quadrupolar line-shape simulation.
"""
#%%
# Coesite is a high-pressure (2-3 GPa) and high-temperature (700°C) polymorph of silicon
# dioxide :math:`\text{SiO}_2`. Coesite has five crystallographic :math:`^{17}\text{O}`
# sites. In the following, use the :math:`^{17}\text{O}` EFG tensor information from
# Grandinetti et. al. [#f2]_
import matplotlib.pyplot as plt
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site

#%%
# **Step 1** Create the sites.

# default unit of isotropic_chemical_shift is ppm and Cq is Hz.
O17_1 = Site(
    isotope="17O", isotropic_chemical_shift=29, quadrupolar={"Cq": 6.05e6, "eta": 0.000}
)
O17_2 = Site(
    isotope="17O", isotropic_chemical_shift=41, quadrupolar={"Cq": 5.43e6, "eta": 0.166}
)
O17_3 = Site(
    isotope="17O", isotropic_chemical_shift=57, quadrupolar={"Cq": 5.45e6, "eta": 0.168}
)
O17_4 = Site(
    isotope="17O", isotropic_chemical_shift=53, quadrupolar={"Cq": 5.52e6, "eta": 0.169}
)
O17_5 = Site(
    isotope="17O", isotropic_chemical_shift=58, quadrupolar={"Cq": 5.16e6, "eta": 0.292}
)

# all five sites.
sites = [O17_1, O17_2, O17_3, O17_4, O17_5]

#%%
# **Step 2** Create the isotopomers from these sites. For optimum performance, we
# create five single-site isotopomers than a five-site single isotopomer.

abundance = [0.83, 1.05, 2.16, 2.05, 1.90]  # abundance of each isotopomer
isotopomers = [Isotopomer(sites=[s], abundance=a) for s, a in zip(sites, abundance)]

#%%
# **Step 3** Create a central transition selective Bloch decay spectrum method.
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum

method = BlochDecayCentralTransitionSpectrum(
    channels=["17O"],
    rotor_frequency=14000,  # in Hz
    spectral_dimensions=[
        {"count": 2046, "spectral_width": 50000}  # spectral_width is in Hz
    ],
)

#%%
# The above method is set up to record the :math:`^{17}\text{O}` resonances at the
# magic angle, spinning at 14 kHz and 9.4 T (default, if the value is not provided)
# external magnetic flux density. The resonances are recorded over 50 kHz spectral
# width using 2046 points.

#%%
# **Step 4** Create the Simulator object and add the method and isotopomer objects.

sim_coesite = Simulator()
sim_coesite.isotopomers += isotopomers  # add isotopomers
sim_coesite.methods += [method]  # add method

#%%
# **Step 5** Simulate the spectrum.

sim_coesite.run()

#%%
# **Step 6** The plot of the simulation.

x, y = sim_coesite.methods[0].simulation.to_list()
plt.figure(figsize=(4, 3))
plt.plot(x, y, color="black", linewidth=1)
plt.xlabel("frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# .. [#f2] Grandinetti, P. J., Baltisberger, J. H., Farnan, I., Stebbins, J. F.,
#       Werner, U. and Pines, A.
#       Solid-State :math:`^{17}\text{O}` Magic-Angle and Dynamic-Angle Spinning NMR
#       Study of the :math:`\text{SiO}_2` Polymorph Coesite, J. Phys. Chem. 1995,
#       **99**, *32*, 12341-12348.
#       `DOI: 10.1021/j100032a045 <https://doi.org/10.1021/j100032a045>`_
