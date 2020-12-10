#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
13C 2D PASS NMR of LHistidine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# Coesite is a high-pressure (2-3 GPa) and high-temperature (700°C) polymorph of silicon
# dioxide :math:`\text{SiO}_2`. Coesite has five crystallographic :math:`^{17}\text{O}`
# sites. The experimental dataset used in this example is published in
# Grandinetti `et. al.` [#f1]_
import numpy as np
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator
from mrsimulator.methods import SSB2D
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params
from lmfit import Minimizer, report_fit
from mrsimulator.utils.collection import single_site_system_generator

# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/687656/files/1H13C_CPPASS_LHistidine.csdf"
pass_data = cp.load(filename)

# For the spectral fitting, we only focus on the real part of the complex dataset.
# The script assumes that the dimension at index 0 is the isotropic dimension.
# Transpose the dataset as required.
pass_data = pass_data.real.T

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in pass_data.dimensions]

# Normalize the spectrum
pass_data /= pass_data.max()

# plot of the dataset.
levels = (np.arange(10) + 0.3) / 15  # contours are drawn at these levels.
ax = plt.subplot(projection="csdm")
cb = ax.contour(pass_data, colors="k", levels=levels, alpha=0.5, linewidths=0.5)
plt.colorbar(cb)
ax.set_xlim(200, 10)
ax.invert_yaxis()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# The fitting model includes the Simulator and the SignalProcessor objects. First
# create the Simulator object.

# Create the guess sites and spin systems.
# default unit of isotropic_chemical_shift is ppm and Cq is Hz.
shifts = [120, 128, 135, 175, 55, 25]  # in ppm
zeta = [-70, -65, -60, -60, -10, -10]  # in  Hz
eta = [0.8, 0.4, 0.9, 0.3, 0.0, 0.0]

spin_systems = single_site_system_generator(
    isotopes="13C",
    isotropic_chemical_shifts=shifts,
    shielding_symmetric={"zeta": zeta, "eta": eta},
    abundance=100 / 6,
)

# Create the DAS method.
# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(pass_data)

# %%
ssb = SSB2D(
    channels=["13C"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=pass_data,  # also add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the das method.
for sys in spin_systems:
    sys.transition_pathways = ssb.get_transition_pathways(sys)

# %%

# Create the Simulator object and add the method and spin system objects.
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [ssb]  # add the method
sim.run()

# %%

# Add Post simulation processing
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along the isotropic dimensions.
        sp.FFT(axis=0),
        apo.Exponential(FWHM="20 Hz"),
        sp.IFFT(axis=0),
        sp.Scale(factor=0.6),
    ]
)
# Apply post simulation operations
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# %%

# The plot of the simulation after signal processing.
ax = plt.subplot(projection="csdm")
ax.contour(processed_data, colors="r", levels=levels, alpha=0.5, linewidths=0.5)
cb = ax.contour(pass_data, colors="k", levels=levels, alpha=0.5, linewidths=0.5)
plt.colorbar(cb)
ax.set_xlim(200, 10)
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# First create the fitting parameters.
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup.
params = make_LMFIT_params(sim, processor)
print(params.pretty_print())

# %%
# Run the minimization using LMFIT
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor))
result = minner.minimize()
report_fit(result)


# %%
# Simulate the spectrum corresponding to the optimum parameters
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# %%
# Plot the spectrum
ax = plt.subplot(projection="csdm")
ax.contour(processed_data, colors="r", levels=levels, alpha=0.5, linewidths=0.5)
cb = ax.contour(pass_data, colors="k", levels=levels, alpha=0.5, linewidths=0.5)
plt.colorbar(cb)
ax.set_xlim(200, 10)
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Grandinetti, P. J., Baltisberger, J. H., Farnan, I., Stebbins, J. F.,
#       Werner, U. and Pines, A.
#       Solid-State :math:`^{17}\text{O}` Magic-Angle and Dynamic-Angle Spinning NMR
#       Study of the :math:`\text{SiO}_2` Polymorph Coesite, J. Phys. Chem. 1995,
#       **99**, *32*, 12341-12348.
#       `DOI: 10.1021/j100032a045 <https://doi.org/10.1021/j100032a045>`_
