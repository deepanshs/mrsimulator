# -*- coding: utf-8 -*-
"""Lineshape Test."""
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.affine as aft
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.methods import Method2D
from mrsimulator.methods import ThreeQ_VAS


def test_MQMAS():
    spin_system = SpinSystem(
        sites=[
            Site(
                isotope="87Rb",
                isotropic_chemical_shift=-9,
                shielding_symmetric={"zeta": 100, "eta": 0},
                quadrupolar={"Cq": 3.5e6, "eta": 0.36, "beta": 70 / 180 * np.pi},
            )
        ]
    )

    method = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        spectral_dimensions=[
            {
                "count": 128,
                "spectral_width": 20000,
                "events": [{"transition_query": {"P": [-3], "D": [0]}}],
            },
            {"count": 128, "spectral_width": 20000},
        ],
    )

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.config.integration_volume = "hemisphere"
    sim.run()

    # process
    k = 21 / 27  # shear factor
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=1),
            aft.Shear(factor=-k, dim_index=1, parallel=0),
            aft.Scale(factor=1 + k, dim_index=1),
            sp.FFT(dim_index=1),
        ]
    )
    processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

    # Since there is a single site, after the shear and scaling transformations, there
    # should be a single perak along the isotropic dimension at index 70.
    # The isotropic coordinate of this peak is given by
    # w_iso = (17.8)*iso_shift + 1e6/8 * (vq/v0)^2 * (eta^2 / 3 + 1)
    # ref: D. Massiot et al. / Solid State Nuclear Magnetic Resonance 6 (1996) 73-83
    iso_slice = processed_data[40, :]
    assert np.argmax(iso_slice.dependent_variables[0].components[0]) == 70

    # calculate the isotropic coordinate
    spin = method.channels[0].spin
    w0 = method.channels[0].gyromagnetic_ratio * 9.4 * 1e6
    wq = 3 * 3.5e6 / (2 * spin * (2 * spin - 1))
    w_iso = -9 * 17 / 8 + 1e6 / 8 * (wq / w0) ** 2 * ((0.36 ** 2) / 3 + 1)

    # the coordinate from spectrum
    w_iso_spectrum = processed_data.dimensions[1].coordinates[70].value
    np.testing.assert_almost_equal(w_iso, w_iso_spectrum, decimal=2)

    # The projection onto the  MAS dimension should be the 1D block decay central
    # transition spectrum
    mas_slice = processed_data.sum(axis=1).dependent_variables[0].components[0]

    # MAS spectrum
    method = BlochDecayCentralTransitionSpectrum(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=1e9,
        spectral_dimensions=[{"count": 128, "spectral_width": 20000}],
    )

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.config.integration_volume = "hemisphere"
    sim.run()

    data = sim.methods[0].simulation.dependent_variables[0].components[0]
    assert np.allclose(data / data.max(), mas_slice / mas_slice.max())


def test_ThreeQ_VAS_spin_3halves():
    spin_system = SpinSystem(
        sites=[
            Site(
                isotope="87Rb",
                isotropic_chemical_shift=-9,
                shielding_symmetric={"zeta": 100, "eta": 0},
                quadrupolar={"Cq": 3.5e6, "eta": 0.36, "beta": 70 / 180 * np.pi},
            )
        ]
    )

    method = ThreeQ_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 20000},
            {"count": 512, "spectral_width": 20000},
        ],
    )
    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.config.integration_volume = "hemisphere"
    sim.run()

    data = sim.methods[0].simulation
    dat = data.dependent_variables[0].components[0]
    index = np.where(dat == dat.max())

    # The isotropic coordinate of this peak is given by
    # v_iso = (17/8)*iso_shift + 1e6/8 * (vq/v0)^2 * (eta^2 / 3 + 1)
    # ref: D. Massiot et al. / Solid State Nuclear Magnetic Resonance 6 (1996) 73-83
    spin = method.channels[0].spin
    v0 = method.channels[0].gyromagnetic_ratio * 9.4 * 1e6
    vq = (3 * 3.5e6) / (2 * spin * (2 * spin - 1))
    v_iso = -9 * 17 / 8 + 1e6 / 8 * ((vq / v0) ** 2) * ((0.36 ** 2) / 3 + 1)

    # the coordinate from spectrum along the iso dimension must be equal to v_iso
    v_iso_spectrum = data.dimensions[1].coordinates[index[0]].value
    np.testing.assert_almost_equal(v_iso, v_iso_spectrum, decimal=2)

    # The projection onto the  MAS dimension should be the 1D block decay central
    # transition spectrum
    mas_slice = data.sum(axis=1).dependent_variables[0].components[0]

    # MAS spectrum
    method = BlochDecayCentralTransitionSpectrum(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=1e9,
        spectral_dimensions=[{"count": 512, "spectral_width": 20000}],
    )

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.config.integration_volume = "hemisphere"
    sim.run()

    data = sim.methods[0].simulation.dependent_variables[0].components[0]
    assert np.allclose(data / data.max(), mas_slice / mas_slice.max())


def test_MQMAS_spin_5halves():
    spin_system = SpinSystem(
        sites=[
            Site(
                isotope="27Al",
                isotropic_chemical_shift=64.5,  # in ppm
                quadrupolar={"Cq": 3.22e6, "eta": 0.66},  # Cq is in Hz
            )
        ]
    )

    method = ThreeQ_VAS(
        channels=["27Al"],
        magnetic_flux_density=7,
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 5000, "reference_offset": -3e3},
            {"count": 512, "spectral_width": 10000, "reference_offset": 4e3},
        ],
    )

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.run()

    data = sim.methods[0].simulation
    dat = data.dependent_variables[0].components[0]
    index = np.where(dat == dat.max())

    # The isotropic coordinate of this peak is given by
    # v_iso = -(17/31)*iso_shift + 8e6/93 * (vq/v0)^2 * (eta^2 / 3 + 1)
    # ref: D. Massiot et al. / Solid State Nuclear Magnetic Resonance 6 (1996) 73-83
    spin = method.channels[0].spin
    v0 = method.channels[0].gyromagnetic_ratio * 7 * 1e6
    vq = 3 * 3.22e6 / (2 * spin * (2 * spin - 1))
    v_iso = -(17 / 31) * 64.5 - (8e6 / 93) * (vq / v0) ** 2 * ((0.66 ** 2) / 3 + 1)

    # the coordinate from spectrum along the iso dimension must be equal to v_iso
    v_iso_spectrum = data.dimensions[1].coordinates[index[0]].value
    np.testing.assert_almost_equal(v_iso, v_iso_spectrum, decimal=2)

    # The projection onto the  MAS dimension should be the 1D block decay central
    # transition spectrum
    mas_slice = data.sum(axis=1).dependent_variables[0].components[0]

    # MAS spectrum
    method = BlochDecayCentralTransitionSpectrum(
        channels=["27Al"],
        magnetic_flux_density=7,
        rotor_frequency=1e9,
        spectral_dimensions=[
            {"count": 512, "spectral_width": 10000, "reference_offset": 4e3}
        ],
    )

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [method]
    sim.config.integration_volume = "hemisphere"
    sim.run()

    data = sim.methods[0].simulation.dependent_variables[0].components[0]
    assert np.allclose(data / data.max(), mas_slice / mas_slice.max())
