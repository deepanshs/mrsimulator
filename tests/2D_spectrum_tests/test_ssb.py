# -*- coding: utf-8 -*-
"""Lineshape Test."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.methods import Method2D
from mrsimulator.methods import SSB2D


def SSB2D_setup(ist, vr, method_type):
    sites = [
        Site(
            isotope=ist,
            isotropic_chemical_shift=29,
            shielding_symmetric={"zeta": -70, "eta": 0.000},
        ),
        Site(
            isotope=ist,
            isotropic_chemical_shift=44,
            shielding_symmetric={"zeta": -96, "eta": 0.166},
        ),
        Site(
            isotope=ist,
            isotropic_chemical_shift=57,
            shielding_symmetric={"zeta": -120, "eta": 0.168},
        ),
    ]
    spin_systems = [SpinSystem(sites=[s]) for s in sites]

    B0 = 11.7
    if method_type == "PASS":
        method = SSB2D(
            channels=[ist],
            magnetic_flux_density=B0,  # in T
            rotor_frequency=vr,
            spectral_dimensions=[
                {
                    "count": 32,
                    "spectral_width": 32 * vr,  # in Hz
                    "label": "Anisotropic dimension",
                },
                # The last spectral dimension block is the direct-dimension
                {
                    "count": 2048,
                    "spectral_width": 2e4,  # in Hz
                    "reference_offset": 5e3,  # in Hz
                    "label": "High speed MAS dimension",
                },
            ],
        )
    else:
        method = Method2D(
            channels=[ist],
            magnetic_flux_density=B0,  # in T
            spectral_dimensions=[
                {
                    "count": 64,
                    "spectral_width": 8e4,  # in Hz
                    "label": "Anisotropic dimension",
                    "events": [{"rotor_angle": 90 * 3.14159 / 180}],
                },
                # The last spectral dimension block is the direct-dimension
                {
                    "count": 2048,
                    "spectral_width": 2e4,  # in Hz
                    "reference_offset": 5e3,  # in Hz
                    "label": "High speed MAS dimension",
                },
            ],
            affine_matrix=[1, -1, 0, 1],
        )
    sim = Simulator()
    sim.spin_systems = spin_systems  # add spin systems
    sim.methods = [method]  # add the method.
    sim.run()

    data_ssb = sim.methods[0].simulation
    dim_ssb = data_ssb.x[0].coordinates.value

    if method_type == "PASS":
        bloch = BlochDecaySpectrum(
            channels=[ist],
            magnetic_flux_density=B0,  # in T
            rotor_frequency=vr,  # in Hz
            spectral_dimensions=[
                {
                    "count": 32,
                    "spectral_width": 32 * vr,  # in Hz
                    "reference_offset": 0,  # in Hz
                    "label": "MAS dimension",
                },
            ],
        )
    else:
        bloch = BlochDecaySpectrum(
            channels=[ist],
            magnetic_flux_density=B0,  # in T
            rotor_frequency=vr,  # in Hz
            rotor_angle=90 * 3.14159 / 180,
            spectral_dimensions=[
                {
                    "count": 64,
                    "spectral_width": 8e4,  # in Hz
                    "reference_offset": 0,  # in Hz
                    "label": "MAS dimension",
                },
            ],
        )

    for i in range(3):
        iso = spin_systems[i].sites[0].isotropic_chemical_shift
        sys = spin_systems[i].copy()
        sys.sites[0].isotropic_chemical_shift = 0
        sim2 = Simulator()
        sim2.spin_systems = [sys]  # add spin systems
        sim2.methods = [bloch]  # add the method.
        sim2.run()

        index = np.where(dim_ssb < iso)[0][-1]

        one_d_section = data_ssb.y[0].components[0][:, index]
        one_d_section /= one_d_section.max()

        one_d_sim = sim2.methods[0].simulation.y[0].components[0]
        one_d_sim /= one_d_sim.max()

        np.testing.assert_almost_equal(one_d_section, one_d_sim, decimal=6)


def test_01():
    SSB2D_setup("29Si", 1500, "PASS")


def test_02():
    SSB2D_setup("31P", 2500, "PASS")


def test_03():
    SSB2D_setup("29Si", 1e9, "MAF")


def test_04():
    SSB2D_setup("31P", 1e9, "MAF")
