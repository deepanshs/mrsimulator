#!/usr/bin/env python
"""
²H (I=1) 2D sideband-sideband correlation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

²H (I=1) 2D NMR CSA-Quad 1st order sideband correlation spectrum simulation.
"""
# %%
# Sideband-sideband NMR correlation simulation of crystalline solid as
# reported by Aleksis and Pell [#f1]_.
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent
from mrsimulator.simulator.sampling_scheme import zcw_averaging

# sphinx_gallery_thumbnail_number = 1

# %%
# Generate the site and spin system objects.
alpha = [0, 30, 60, 90, 0, 0, 0, 0, 90, 90, 90, 90]
beta = [0, 0, 0, 0, 0, 30, 60, 90, 90, 90, 90, 90]
gamma = [0, 0, 0, 0, 90, 90, 90, 90, 0, 30, 60, 90]

spin_systems = [
    SpinSystem(
        name=f"$\\alpha={al}, \\beta={be}, \\gamma={ga}$",
        sites=[
            Site(
                isotope="2H",
                isotropic_chemical_shift=0,  # in ppm
                shielding_symmetric=SymmetricTensor(
                    zeta=150,  # in ppm
                    eta=0.7,
                    alpha=al * 3.14159 / 180,  # in rads
                    beta=be * 3.14159 / 180,  # in rads
                    gamma=ga * 3.14159 / 180,  # in rads
                ),
                quadrupolar=SymmetricTensor(Cq=50e3, eta=0.9),  # Cq in Hz
            )
        ],
    )
    for al, be, ga in zip(alpha, beta, gamma)
]

# %%
# Create a sideband-sideband correlation method
rotor_frequency = 2000

sideband_2d = Method(
    name="2D sideband correlation",
    channels=["2H"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=rotor_frequency,  # in Hz
    spectral_dimensions=[
        SpectralDimension(
            count=16,
            spectral_width=rotor_frequency * 16,  # in Hz
            label="Paramagnetic shift",
            events=[
                SpectralEvent(
                    transition_queries=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Shielding1_2"],
                ),
                MixingEvent(),
            ],
        ),
        SpectralDimension(
            count=50,
            spectral_width=50 * rotor_frequency,  # in Hz
            label="Quadrupolar + iso frequency",
            events=[
                SpectralEvent(
                    transition_queries=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Quad1_2", "Shielding1_0"],
                )
            ],
        ),
    ],
)

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator(spin_systems=spin_systems, methods=[sideband_2d])
sim.config.decompose_spectrum = "spin_system"  # simulate spectra per spin system
sim.config.number_of_sidebands = 56

# custom sampling scheme
sim.config.custom_sampling = zcw_averaging(
    M=11, integration_volume="hemisphere", triangle_mesh=False
)

sim.run()


# %%
dataset = sim.methods[0].simulation.real
_ = [dim.to("kHz", "nmr_frequency_ratio") for dim in dataset.dimensions]
datasets = dataset.split()

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
cb = ax.imshow(
    datasets[0] / datasets[0].max(),
    aspect="auto",
    cmap="gist_ncar_r",
    interpolation="none",
)
plt.title(None)
plt.colorbar(cb)
plt.tight_layout()
plt.show()

# %%
# Let's plot all the simulated datasets.
fig, ax = plt.subplots(
    3, 4, sharex=True, sharey=True, figsize=(12, 9.5), subplot_kw={"projection": "csdm"}
)

vmax = max(dataset.max())
vmin = min(dataset.min())
for j, datum in enumerate(datasets):
    row, col = j // 4, j % 4
    ax[row, col].imshow(
        datum,
        aspect="auto",
        cmap="gist_ncar_r",
        interpolation="none",
        vmax=vmax,
        vmin=vmin,
    )
    ax[row, col].invert_xaxis()

plt.tight_layout()
plt.show()

# %%
# .. [#f1] Aleksis R and Pell A.J., Separation of quadrupolar and paramagnetic shift
#        interactions in high-resolution nuclear magnetic resonance of spinning
#        powders, J. Chem. Phys. (2021)  **155**, 094202.
#        `DOI: 10.1063/5.0061611 <https://doi.org/10.1063/5.0061611>`_
