"""Lineshape Test."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.spin_system.isotope import Isotope

# default unit of isotropic_chemical_shift is ppm and Cq is Hz.
O17 = Isotope(symbol="17O")
O17_1 = Site(
    isotope="17O",
    isotropic_chemical_shift=29,
    quadrupolar={"Cq": 6.05e6, "eta": 0.000},
)
O17_2 = Site(
    isotope="17O",
    isotropic_chemical_shift=41,
    quadrupolar={"Cq": 5.43e6, "eta": 0.166},
)
O17_3 = Site(
    isotope="17O",
    isotropic_chemical_shift=57,
    quadrupolar={"Cq": 5.45e6, "eta": 0.168},
)
O17_4 = Site(
    isotope="17O",
    isotropic_chemical_shift=53,
    quadrupolar={"Cq": 5.52e6, "eta": 0.169},
)
O17_5 = Site(
    isotope="17O",
    isotropic_chemical_shift=58,
    quadrupolar={"Cq": 5.16e6, "eta": 0.292},
)

# all five sites.
sites = [O17_1, O17_2, O17_3, O17_4, O17_5]
abundance = [1, 1, 2, 2, 2]
spin_systems = [SpinSystem(sites=[s], abundance=a) for s, a in zip(sites, abundance)]


def quad_iso_shift(m, C_q, eta_q, spin, delta_iso_cs, nu_0, nu_ref):
    nu_q = (3 * C_q) / (2 * spin * (2 * spin - 1))
    nu_ref = nu_ref * (-nu_0 / abs(nu_0))  # nu_ref has opposite sign to nu_0
    term1 = 2 * m * delta_iso_cs * (-nu_0 * 1e-6)
    term2 = (2 * m / 30) * ((eta_q**2 / 3) + 1) * (spin * (spin + 1) - 3 * m**2)
    shift = term1 + (nu_q**2 / nu_0) * term2
    # divide shift by nu_ref before returning to convert from Hz to ppm
    return shift


def test_DAS():
    B0 = 11.7
    das = Method(
        channels=["17O"],
        magnetic_flux_density=B0,  # in T
        rotor_frequency=1e12,
        spectral_dimensions=[
            {
                "count": 912,
                "spectral_width": 5e3,  # in Hz
                "reference_offset": 0,  # in Hz
                "origin_offset": O17.B0_to_ref_freq(B0) * 1e6,  # in Hz
                "label": "DAS isotropic dimension",
                "events": [
                    {
                        "fraction": 0.5,
                        "rotor_angle": 37.38 * 3.14159 / 180,
                        "transition_queries": [{"ch1": {"P": [-1], "D": [0]}}],
                    },
                    {
                        "fraction": 0.5,
                        "rotor_angle": 79.19 * 3.14159 / 180,
                        "transition_queries": [{"ch1": {"P": [-1], "D": [0]}}],
                    },
                ],
            },
            # The last spectral dimension block is the direct-dimension
            {
                "count": 2048,
                "spectral_width": 2e4,  # in Hz
                "reference_offset": 0,  # in Hz
                "origin_offset": O17.B0_to_ref_freq(B0) * 1e6,  # in Hz
                "label": "MAS dimension",
                "events": [
                    {
                        "rotor_angle": 54.735 * 3.14159 / 180,
                        "transition_queries": [{"ch1": {"P": [-1], "D": [0]}}],
                    }
                ],
            },
        ],
    )

    sim = Simulator()
    sim.spin_systems = spin_systems  # add spin systems
    sim.methods = [das]  # add the method.
    sim.config.decompose_spectrum = "spin_system"
    sim.run(pack_as_csdm=False)

    dataset_das = sim.methods[0].simulation
    dataset_das_coords_Hz = das.spectral_dimensions[0].coordinates_Hz()
    # dataset_das_coords_ppm = das.spectral_dimensions[0].coordinates_ppm()

    # Bloch decay central transition method
    bloch = BlochDecayCTSpectrum(
        channels=["17O"],
        magnetic_flux_density=B0,  # in T
        rotor_frequency=np.inf,  # in Hz
        rotor_angle=54.735 * 3.14159 / 180,
        spectral_dimensions=[
            {
                "count": 2048,
                "spectral_width": 2e4,  # in Hz
                "reference_offset": 0,  # in Hz
                "label": "MAS dimension",
            },
        ],
    )

    sim = Simulator()
    sim.spin_systems = spin_systems
    sim.methods = [bloch]
    sim.config.decompose_spectrum = "spin_system"
    sim.run(pack_as_csdm=False)

    dataset_bloch = sim.methods[0].simulation

    larmor_freq = das.channels[0].gyromagnetic_ratio * B0 * 1e6
    spin = das.channels[0].spin
    for i, sys in enumerate(spin_systems):
        for site in sys.sites:
            Cq = 0
            eta = 0
            if site.quadrupolar is not None:
                Cq = site.quadrupolar.Cq
                eta = site.quadrupolar.eta
            iso = site.isotropic_chemical_shift
            factor1 = -(3 / 40) * (Cq / larmor_freq) ** 2
            factor2 = (spin * (spin + 1) - 3 / 4) / (spin**2 * (2 * spin - 1) ** 2)
            factor3 = 1 + (eta**2) / 3
            iso_obs = factor1 * factor2 * factor3 * 1e6 + iso
            # iso_obs = quad_iso_shift(
            #     0.5, Cq, eta, spin, iso, larmor_freq, O17.B0_to_ref_freq(B0)
            # )

            # get the index where there is a signal
            id1 = dataset_das[i] / dataset_das[i].max()
            index = np.where(id1 == id1.max())[0]
            iso_spectrum = dataset_das_coords_Hz[index[0]]  # x[1].coords[index[0]]
            # iso_spectrum = dataset_das_coords_ppm[index[0]]  # x[1].coords[index[0]]

            print(iso_obs, iso_spectrum)
            # test for the position of isotropic peaks.
            np.testing.assert_almost_equal(iso_obs, iso_spectrum, decimal=1)

            # test for the spectrum across the isotropic peaks.
            dataset_bloch_i = dataset_bloch[i] / dataset_bloch[i].max()
            assert np.allclose(id1[index[0]], dataset_bloch_i)
