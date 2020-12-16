#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
11B MAS NMR of Lithium orthoborate crystal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a quadrupolar lineshape fitting example for the 11B MAS NMR of
# lithium orthoborate crystal. The dataset was provided by Nathan Barrow.
import numpy as np
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator, Site, SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params
from lmfit import Minimizer, report_fit

# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
mpl.rcParams["grid.linestyle"] = "--"
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/710705/files/11B_lithum_orthoborate.csdf"
experiment = cp.load(filename)

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# Normalize the spectrum
experiment /= experiment.max()

# plot of the dataset.
levels = (np.arange(10) + 0.3) / 15  # contours are drawn at these levels.
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
ax.set_xlim(50, -25)
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.
B11 = Site(
    isotope="11B",
    isotropic_chemical_shift=20.0,  # in ppm
    quadrupolar={"Cq": 2.3e6, "eta": 0.03},  # Cq in Hz
)
spin_systems = [SpinSystem(sites=[B11])]


# %%
# **Method**

# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecayCentralTransitionSpectrum(
    channels=["11B"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=12500,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the das method.
for sys in spin_systems:
    sys.transition_pathways = method.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Lorentzian convolution along both dimensions.
        sp.IFFT(),
        apo.Exponential(FWHM="100 Hz"),
        sp.FFT(),
        sp.Scale(factor=1),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5, linewidth=2, label="Experiment")
ax.plot(processed_data, "r", label="guess spectrum")
ax.set_xlim(50, -25)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = make_LMFIT_params(sim, processor)
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor))
result = minner.minimize()
report_fit(result)


# %%
# The best fit solution
# ---------------------
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot the spectrum
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5, linewidth=2, label="Experiment")
ax.plot(processed_data, "r--", label="Best Fit")
ax.set_xlim(50, -25)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
