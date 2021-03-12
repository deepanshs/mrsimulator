# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.utils import flatten_dict

from .base import ContribSchema

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


SITE_KEYWORDS = {
    "isotropic_chemical_shift": "isotropic",
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
    rotor_frequency = method.spectral_dimensions[0].events[0].rotor_frequency  # Hz
    spectral_width = method.spectral_dimensions[0].spectral_width  # Hz
    rotor_angle = method.spectral_dimensions[0].events[0].rotor_angle * 180 / np.pi

    return {
        "larmorFrequency": f"{larmor_frequency} MHz",
        "spinningFrequency": f"{rotor_frequency} Hz",
        "spectralWidth": f"{spectral_width} Hz",
        "rotorAngle": f"{rotor_angle:.4f} degree",
    }


def mpcontribs_export(sim, project, identifier, exp_dict={}, **kwargs):
    """Generate mpcontribs contribution entries for every site in the Simulator object.

    Arguments
    ---------
        sim: Simulator object from where the site contributions are extracted.
        str project: mpcontribs project name (reqiuired).
        str identifier: mpcontribs identifier (required).
        exp_dict: Additional metadata to use in contribution.
        **kwargs: Optional keyword arguments from mpcontribs ContributionsSchema

    Example
    -------
        >>> contribution_data = mpcontribs_export(sim, 'myproject') # doctest:+SKIP
    """
    contribs = [
        ContribSchema(
            project=project,
            identifier=identifier,
            data={
                "experiment": "experiment goes here",
                "simulation": "simulation goes here",
                "site": {**parse_sites(site)},
                "method": {**parse_method(sim.methods[0])},
            },
            **kwargs,
        ).json()
        for sys in sim.spin_systems
        for site in sys.sites
    ]

    for item in contribs:
        item["data"]["method"] = {**item["data"]["method"], **exp_dict}

    return contribs
