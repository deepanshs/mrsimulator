#!/usr/bin/env python
"""
MCl₂.2D₂O, ²H (I=1) Shifting-d echo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

²H (I=1) 2D NMR CSA-Quad 1st order sideband correlation spectrum simulation.
"""
# %%
# Sideband-sideband NMR correlation simulation of crystalline solid as
# reported by Aleksis and Pell [#f1]_.
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent

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
n_sidebands = 48
rotor_frequency = 2000
bandwidth = rotor_frequency * n_sidebands

shifting_d = Method(
    name="Shifting-d",
    channels=["2H"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=rotor_frequency,  # in Hz
    spectral_dimensions=[
        SpectralDimension(
            count=n_sidebands,
            spectral_width=bandwidth,  # in Hz
            label="Quadrupolar frequency",
            events=[
                SpectralEvent(
                    transition_queries=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Quad1_2"],
                ),
                MixingEvent(query="NoMixing"),
            ],
        ),
        SpectralDimension(
            count=n_sidebands,
            spectral_width=bandwidth,  # in Hz
            label="Paramagnetic shift",
            events=[
                SpectralEvent(
                    transition_queries=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Shielding1_0", "Shielding1_2"],
                )
            ],
        ),
    ],
)

# %%
# Create the Simulator object, add the method and spin system objects, and
# run the simulation.
sim = Simulator(spin_systems=spin_systems, methods=[shifting_d])
sim.config.integration_volume = "hemisphere"
sim.config.decompose_spectrum = "spin_system"  # simulate spectra per spin system
sim.config.number_of_sidebands = n_sidebands
sim.run()

# %%
# Let's plot all the simulated datasets.
dataset = sim.methods[0].simulation
[dim.to("kHz", "nmr_frequency_ratio") for dim in dataset.dimensions]
fig, ax = plt.subplots(
    3, 4, sharex=True, sharey=True, figsize=(12, 9.5), subplot_kw={"projection": "csdm"}
)
for j, datum in enumerate(dataset.split()):
    row, col = j // 4, j % 4
    ax[row, col].imshow(
        (datum.T / datum.max()).real,
        aspect="auto",
        cmap="gist_ncar_r",
        interpolation="none",
    )
    ax[row, col].set_ylim(-14.5, 14.5)
    ax[row, col].invert_xaxis()
    ax[row, col].invert_yaxis()

plt.tight_layout()
plt.show()

# %%
# .. [#f1] Aleksis R and Pell A.J., Separation of quadrupolar and paramagnetic shift
#        interactions in high-resolution nuclear magnetic resonance of spinning
#        powders, J. Chem. Phys. (2021)  **155**, 094202.
#        `DOI: 10.1063/5.0061611 <https://doi.org/10.1063/5.0061611>`_
