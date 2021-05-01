# -*- coding: utf-8 -*-
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.contribs import mpcontribs_export
from mrsimulator.contribs import parse_method
from mrsimulator.contribs import parse_sites
from mrsimulator.methods import BlochDecaySpectrum


def test01():
    site = Site(isotope="1H", shielding_symmetric={"zeta": 10, "eta": 0.1})
    method = BlochDecaySpectrum(
        channels=["1H"], magnetic_flux_density=9.4, rotor_frequency="15000"
    )

    output_site = parse_sites(site)
    output_method = parse_method(method)

    assert output_site == {
        "isotope": "1H",
        "ChemicalShift": {"isotropic": "0.0 ppm", "zeta": "-10.0 ppm", "eta": 0.1},
    }

    omega_0 = abs(method.channels[0].gyromagnetic_ratio * 9.4)
    assert output_method == {
        "larmorFrequency": f"{omega_0} MHz",
        "spinningFrequency": "15000.0 Hz",
        "spectralWidth": "25000.0 Hz",
        "rotorAngle": "54.7356 degree",
    }


def test02():
    site = Site(isotope="27Al", quadrupolar={"Cq": 10e6, "eta": 0.4})
    method = BlochDecaySpectrum(channels=["27Al"], magnetic_flux_density=11.7)

    output_site = parse_sites(site)
    output_method = parse_method(method)

    assert output_site == {
        "isotope": "27Al",
        "ChemicalShift": {"isotropic": "0.0 ppm"},
        "Quadrupolar": {"Cq": "10.0 MHz", "eta": 0.4},
    }

    omega_0 = abs(method.channels[0].gyromagnetic_ratio * 11.7)
    assert output_method == {
        "larmorFrequency": f"{omega_0} MHz",
        "spinningFrequency": "0.0 Hz",
        "spectralWidth": "25000.0 Hz",
        "rotorAngle": "54.7356 degree",
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
            "site": {
                "isotope": "27Al",
                "ChemicalShift": {"isotropic": "0.0 ppm"},
                "Quadrupolar": {"Cq": "10.0 MHz", "eta": 0.4},
            },
            "method": {
                "larmorFrequency": f"{omega_0} MHz",
                "spinningFrequency": "0.0 Hz",
                "spectralWidth": "25000.0 Hz",
                "rotorAngle": "54.7356 degree",
                "blah": "blah",
            },
        },
        "project": "test",
        "formula": "ABX",
        "identifier": "blah-blah",
    }

    output = mpcontribs_export(
        sim,
        project="test",
        identifier="blah-blah",
        exp_dict={"blah": "blah"},
        formula="ABX",
    )
    assert output == [card]

    sim.spin_systems = [SpinSystem(sites=[site, site, site])]
    output = mpcontribs_export(
        sim,
        project="test",
        identifier="mp-5733",
        exp_dict={"blah": "blah"},
        formula="ABX",
    )
    card["identifier"] = "mp-5733"
    assert output == [card, card, card]
