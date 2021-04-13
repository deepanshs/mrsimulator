#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
119Sn MAS NMR of SnO
^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a spinning sideband manifold fitting example for the 119Sn MAS NMR
# of SnO. The dataset was acquired and shared by Altenhof `et. al.` [#f1]_.
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo

from mrsimulator import Simulator, SpinSystem, Site, Coupling
from mrsimulator.methods import BlochDecaySpectrum
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
filename = "https://sandbox.zenodo.org/record/743242/files/119Sn_SnO.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 0.002697238

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# Normalize the spectrum
max_amp = experiment.max()
experiment /= max_amp
sigma /= max_amp

# plot of the dataset.
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
ax.set_xlim(-1200, 600)
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems. There are two spin systems present in this
# example,
# - 1) an uncoupled :math:`^{119}\text{Sn}` and
# - 2) a coupled :math:`^{119}\text{Sn}`-:math:`^{117}\text{Sn}` spin systems.
sn119 = Site(
    isotope="119Sn",
    isotropic_chemical_shift=-210,
    shielding_symmetric={"zeta": 700, "eta": 0.1},
)
sn117 = Site(
    isotope="117Sn",
    isotropic_chemical_shift=0,
)
j_sn = Coupling(
    site_index=[0, 1],
    isotropic_j=8150.0,
)

sn117_abundance = 7.68  # in %
spin_systems = [
    # uncoupled spin system
    SpinSystem(sites=[sn119], abundance=100 - sn117_abundance),
    # coupled spin systems
    SpinSystem(sites=[sn119, sn117], couplings=[j_sn], abundance=sn117_abundance),
]

# %%
# **Method**

# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecaySpectrum(
    channels=["119Sn"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=10000,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
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
        # Lorentzian convolution.
        sp.IFFT(),
        apo.Exponential(FWHM="1500 Hz"),
        sp.FFT(),
        sp.Scale(factor=20),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.5, linewidth=2.5, label="guess spectrum")
ax.set_xlim(-1200, 600)
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

# Remove the abundance parameters from params. Since the measurement detects 119Sn, we
# also remove the isotropic chemical shift parameter of 117Sn site from params. The
# 117Sn is the site at index 1 of the spin system at index 1.
params.pop("sys_0_abundance")
params.pop("sys_1_abundance")
params.pop("sys_1_site_1_isotropic_chemical_shift")

# Since the 119Sn site is shared between the two spin systems, we add constraints to the
# 119Sn site parameters from the spin system at index 1 to be the same as 119Sn site
# parameters from the spin system at index 0.
lst = [
    "isotropic_chemical_shift",
    "shielding_symmetric_zeta",
    "shielding_symmetric_eta",
]
for item in lst:
    params[f"sys_1_site_0_{item}"].expr = f"sys_0_site_0_{item}"

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
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.5, linewidth=2.5, label="Best Fit")
ax.set_xlim(-1200, 600)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# %%
#
# .. [#f1] Altenhof A. R., Jaroszewicz M. J., Lindquist A. W., Foster L. D. D.,
#       Veinberg S. L., and Schurko R. W. Practical Aspects of Recording Ultra-Wideline
#       NMR Patterns under Magic-Angle Spinning Conditions.
#       J. Phys. Chem. C. 2020, **124**, 27, 14730–14744
#       `DOI: 10.1021/acs.jpcc.0c04510  <https://doi.org/10.1021/acs.jpcc.0c04510>`_
