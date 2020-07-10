#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Post Simulation Signal Processing.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. sectionauthor:: Maxwell C. Venetos <maxvenetos@gmail.com>
"""
# %%
# After running a simulation, we often want to process the resulting spectrum.
# For example, we may want to scale the intensities to match the experiment or
# apodize the signal to simulate line broadening. The following example will
# demonstrate the use of the `SignalProcessor` class to apply various operations
# to simulation data.
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum

# global plot configuration
font = {"weight": "light", "size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 4

# %%
# We will create a hypothetical two-site Si simulation to illustrate post-simulation
# signal processing. We will begin by processing the entire spectrum and follow up by
# decomposing the spectrum and processing each signal independently.
#
# **Step 1** Create the sites and spin system

site1 = {
    "isotope": "29Si",
    "isotropic_chemical_shift": "-75.0 ppm",
    "shielding_symmetric": {"zeta": "-60 ppm", "eta": 0.6},
}

site2 = {
    "isotope": "29Si",
    "isotropic_chemical_shift": "-80.0 ppm",
    "shielding_symmetric": {"zeta": "-70 ppm", "eta": 0.5},
}

spin_system_1 = {
    "name": "site A",
    "description": "A test 29Si site",
    "sites": [site1],  # from the above code
    "abundance": "100%",
}

spin_system_2 = {
    "name": "site B",
    "description": "A test 29Si site",
    "sites": [site2],  # from the above code
    "abundance": "100%",
}

system_object_1 = SpinSystem.parse_dict_with_units(spin_system_1)
system_object_2 = SpinSystem.parse_dict_with_units(spin_system_2)

# %%
# **Step 2** Create the simulation and add the spin system objects.
sim = Simulator()
sim.spin_systems += [system_object_1, system_object_2]

# %%
# **Step 3** Create a Bloch decay spectrum method.
method_dict = {
    "channels": ["29Si"],
    "magnetic_flux_density": "7.1 T",
    "rotor_angle": "54.735 deg",
    "rotor_frequency": "780 Hz",
    "spectral_dimensions": [
        {"count": 2048, "spectral_width": "25 kHz", "reference_offset": "-5 kHz"}
    ],
}
method_object = BlochDecaySpectrum.parse_dict_with_units(method_dict)
sim.methods += [method_object]

# %%
# The above method is set up to record the :math:`^{29}\text{Si}` resonances at the
# magic angle, spinning at 780 Hz and 7.1 T external magnetic flux density. The
# resonances are recorded over 25 kHz spectral width using 2048 points.
#
# **Step 4** Simulate the spectra.
sim.run()

# %%
# **Step 5** Plot the spectrum.
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# **Step 6** Create the post simulation operation by first creating a list of
# processing operations

# list of processing operations
op_list1 = [
    sp.IFFT(dim_indx=0),
    apo.Gaussian(sigma=100, dim_indx=0, dep_var_indx=0),
    sp.FFT(dim_indx=0),
]

# %%
# and then create the :class:`~mrsimulator.signal_processing.SignalProcessor` class
# object as follows
post_sim = sp.SignalProcessor(data=sim.methods[0].simulation, operations=op_list1)

# %%
# The above signal processing procedure will apply the list of operations
# sequentially to the 0-index dependent variable to result in Gaussian
# broadening of the spectrum. The operation procedure will first perform
# a Fourier Transform to convert the data in frequency space to time space.
# Next, a Gaussian exponential with a broadening factor of 100 seconds is
# multipled to the data and then finally another Fourier transform is applied
# to the data to convert the data back to frequency space.
#
#
# **Step 7**  Apply the signal processing.
processed_data = post_sim.apply_operations()

# %%
# **Step 8** Plot the processed spectrum.
ax = plt.subplot(projection="csdm")
ax.plot(processed_data, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# The above code resulted in the same processing to be applied
# to both signals because in the simulation the signals were not
# seperated. In order to apply different processes to each signal,
# we must set the simulation config to decompose the spectrum.
# Steps 1-3 will be the same and we will start at step 4.
#
# **Step 4** Decompose spectrum and run simulation.
sim.config.decompose_spectrum = "spin_system"
sim.run()

# %%
# **Step 5** Plot spectrum.
x, y0, y1 = sim.methods[0].simulation.to_list()

plt.plot(x, y0, x, y1)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

# %%
# **Step 6** Create Post Simulation.
op_list2 = [
    sp.IFFT(dim_indx=0),
    apo.Gaussian(sigma=100, dim_indx=0, dep_var_indx=0),
    apo.Exponential(Lambda=200, dim_indx=0, dep_var_indx=1),
    sp.FFT(dim_indx=0),
]

post_sim = sp.SignalProcessor(data=sim.methods[0].simulation, operations=op_list2)

# %%
# The above signal processing procedure will apply Gaussian
# broadening to dependent variable index 0 (-75 ppm Si site)
# and apply Lorentzian broadening to dependent variable
# index 1 (-80 ppm Si site).

# %%
# **Step 7** Apply the singla processing.
processed_data = post_sim.apply_operations()

# %%
# **Step 8** Plot the processed spectrum
x, y0, y1 = processed_data.to_list()

plt.plot(x, y0, x, y1)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()
