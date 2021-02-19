# -*- coding: utf-8 -*-
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.utils.contribs import contribs_data
from mrsimulator.utils.contribs import parse_method
from mrsimulator.utils.contribs import parse_sites


def test01():
    site = Site(isotope="1H", shielding_symmetric={"zeta": 10, "eta": 0.1})
    method = BlochDecaySpectrum(
        channels=["1H"], magnetic_flux_density=9.4, rotor_frequency="15000"
    )

    output_site = parse_sites(site)
    output_method = parse_method(method)

    assert output_site == {
        "isotope": "1H",
        "ChemicalShift": {"Isotropic": "0 ppm", "zeta": "-10.0 ppm", "eta": 0.1},
    }

    omega_0 = abs(method.channels[0].gyromagnetic_ratio * 9.4)
    assert output_method == {
        "LarmorFrequency": f"{omega_0} MHz",
        "SpinningFrequency": "15000.0 Hz",
        "SpectralWidth": "25000.0 Hz",
        "RotorAngle": "54.7356 degree",
    }


def test02():
    site = Site(isotope="27Al", quadrupolar={"Cq": 10e6, "eta": 0.4})
    method = BlochDecaySpectrum(channels=["27Al"], magnetic_flux_density=11.7)

    output_site = parse_sites(site)
    output_method = parse_method(method)

    assert output_site == {
        "isotope": "27Al",
        "ChemicalShift": {"Isotropic": "0 ppm"},
        "Quadrupolar": {"Cq": "10.0 MHz", "eta": 0.4},
    }

    omega_0 = abs(method.channels[0].gyromagnetic_ratio * 11.7)
    assert output_method == {
        "LarmorFrequency": f"{omega_0} MHz",
        "SpinningFrequency": "0.0 Hz",
        "SpectralWidth": "25000.0 Hz",
        "RotorAngle": "54.7356 degree",
    }


def test_contrib_card():
    site = Site(isotope="27Al", quadrupolar={"Cq": 10e6, "eta": 0.4})
    method = BlochDecaySpectrum(channels=["27Al"], magnetic_flux_density=11.7)

    sim = Simulator()
    sim.spin_systems = [SpinSystem(sites=[site])]
    sim.methods = [method]

    omega_0 = abs(method.channels[0].gyromagnetic_ratio * 11.7)
    card = {
        "data": {
            "experiment": "experiment goes here",
            "simulation": "simulation goes here",
            "site": {
                "isotope": "27Al",
                "ChemicalShift": {"Isotropic": "0 ppm"},
                "Quadrupolar": {"Cq": "10.0 MHz", "eta": 0.4},
            },
            "method": {
                "LarmorFrequency": f"{omega_0} MHz",
                "SpinningFrequency": "0.0 Hz",
                "SpectralWidth": "25000.0 Hz",
                "RotorAngle": "54.7356 degree",
            },
        },
        "project": "test",
    }

    output = contribs_data(sim, "test")
    assert output == [card]

    sim.spin_systems = [SpinSystem(sites=[site, site, site])]
    output = contribs_data(sim, "test", identifier="mp-5733")
    card["identifier"] = "mp-5733"
    assert output == [card, card, card]
