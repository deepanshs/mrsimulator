# -*- coding: utf-8 -*-
import numpy as np

from . import flatten_dict

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


SITE_KEYWORDS = {
    "isotropic_chemical_shift": "Isotropic",
    "shielding_symmetric.zeta": "zeta",
    "shielding_symmetric.eta": "eta",
    "quadrupolar.Cq": "Cq",
    "quadrupolar.eta": "eta",
}


def parse_sites(site):
    if site.shielding_symmetric is not None:
        site.shielding_symmetric.zeta *= -1
    site_ = flatten_dict(site.json())
    dict_ = {"ChemicalShift": {}, "Quadrupolar": {}}
    for k, v in site_.items():
        new_key = SITE_KEYWORDS[k] if k in SITE_KEYWORDS else k
        if "shielding_symmetric" in k:
            dict_["ChemicalShift"][new_key] = site_[k]
        elif "quadrupolar" in k:
            dict_["Quadrupolar"][new_key] = site_[k]
        elif "isotropic_chemical_shift" in k:
            dict_["ChemicalShift"][new_key] = v
        else:
            dict_[new_key] = v

    for item in ["ChemicalShift", "Quadrupolar"]:
        if dict_[item] == {}:
            dict_.pop(item)
    return dict_


def parse_method(method):
    gamma = method.channels[0].gyromagnetic_ratio
    B0 = method.spectral_dimensions[0].events[0].magnetic_flux_density
    larmor_frequency = abs(gamma * B0)  # MHz
    nu_r = method.spectral_dimensions[0].events[0].rotor_frequency  # Hz
    nu_delta = method.spectral_dimensions[0].spectral_width  # Hz
    rotor_angle = method.spectral_dimensions[0].events[0].rotor_angle * 180 / np.pi

    return {
        "LarmorFrequency": f"{larmor_frequency} MHz",
        "SpinningFrequency": f"{nu_r} Hz",
        "SpectralWidth": f"{nu_delta} Hz",
        "RotorAngle": f"{rotor_angle:.4f} degree",
    }


def contribs_data(sim, project, composition=None, identifier=None, exp_dict={}):
    """Generate mpcontribs cards for every site in the Simulator object.

    Arguments
    ---------
        sim: Simulator object from where the sites are extracted.
        project: mpcontribs project name.
        composition: Chemical composition for the sample (optional).
        identifier: The mp-id of the sample (optional).
        exp_dict: Additional metadata to use in contribs card.

    Example
    -------

        >>> contribution_data = contribs_data(sim, 'myproject') # doctest:+SKIP
    """
    contrib = []
    for sys in sim.spin_systems:
        for site in sys.sites:
            data = {
                "experiment": "experiment goes here",
                "simulation": "simulation goes here",
                "site": {**parse_sites(site)},
                "method": {**parse_method(sim.methods[0]), **exp_dict},
            }

            card = {
                "data": data,
                "project": project,
                "composition": composition,
                "identifier": identifier,
            }
            if composition is None:
                card.pop("composition")
            if identifier is None:
                card.pop("identifier")
            contrib.append(card)

    return contrib
