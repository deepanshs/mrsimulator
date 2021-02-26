#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
17O 2D DAS NMR of Coesite
^^^^^^^^^^^^^^^^^^^^^^^^^
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
from mrsimulator.methods import Method2D
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.collection import single_site_system_generator
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params
from lmfit import Minimizer, report_fit


# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
mpl.rcParams["lines.linewidth"] = 0.5
mpl.rcParams["grid.linestyle"] = "--"
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/687656/files/DASCoesite.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 0.1888026

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# Normalize the spectrum
max_amp = experiment.max()
experiment /= experiment.max()
sigma /= max_amp

# plot of the dataset.
levels = (np.arange(10) + 0.3) / 15  # contours are drawn at these levels.
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", levels=levels, alpha=0.75)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.

shifts = [29, 41, 57, 53, 58]  # in ppm
Cq = [6.1e6, 5.4e6, 5.5e6, 5.5e6, 5.1e6]  # in  Hz
eta = [0.1, 0.2, 0.1, 0.1, 0.3]
abundance = [1, 1, 2, 2, 2]

spin_systems = single_site_system_generator(
    isotopes="17O",
    isotropic_chemical_shifts=shifts,
    quadrupolar={"Cq": Cq, "eta": eta},
    abundance=abundance,
)

# %%
# **Method**

# Create the DAS method.
# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

das = Method2D(
    channels=["17O"],
    magnetic_flux_density=11.7,  # in T
    spectral_dimensions=[
        {
            **spectral_dims[0],
            "events": [
                {
                    "fraction": 0.5,
                    "rotor_angle": 37.38 * 3.14159 / 180,
                    "transition_query": {"P": [-1], "D": [0]},
                },
                {
                    "fraction": 0.5,
                    "rotor_angle": 79.19 * 3.14159 / 180,
                    "transition_query": {"P": [-1], "D": [0]},
                },
            ],
        },
        # The last spectral dimension block is the direct-dimension
        {
            **spectral_dims[1],
            "events": [
                {
                    "rotor_angle": 54.735 * 3.14159 / 180,
                    "transition_query": {"P": [-1], "D": [0]},
                }
            ],
        },
    ],
    experiment=experiment,  # also add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the das method.
for sys in spin_systems:
    sys.transition_pathways = das.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [das]  # add the method
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        apo.Gaussian(FWHM="0.15 kHz", dim_index=0),
        apo.Gaussian(FWHM="0.15 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
        sp.Scale(factor=1 / 8),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", levels=levels, alpha=0.75)
ax.contour(processed_data, colors="r", linestyles="--", levels=levels, alpha=0.75)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.grid()
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
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
report_fit(result)


# %%
# The best fit solution
# ---------------------
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot the spectrum
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", levels=levels, alpha=0.75)
ax.contour(processed_data, colors="r", linestyles="--", levels=levels, alpha=0.75)
ax.invert_xaxis()
ax.set_ylim(30, -30)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] Grandinetti, P. J., Baltisberger, J. H., Farnan, I., Stebbins, J. F.,
#       Werner, U. and Pines, A.
#       Solid-State :math:`^{17}\text{O}` Magic-Angle and Dynamic-Angle Spinning NMR
#       Study of the :math:`\text{SiO}_2` Polymorph Coesite, J. Phys. Chem. 1995,
#       **99**, *32*, 12341-12348.
#       `DOI: 10.1021/j100032a045 <https://doi.org/10.1021/j100032a045>`_
