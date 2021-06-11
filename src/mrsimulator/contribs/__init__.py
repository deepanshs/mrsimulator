# -*- coding: utf-8 -*-
import json

import numpy as np
from monty.io import zopen
from mpcontribs.client import Attachment
from mrsimulator import dict
from mrsimulator import Simulator
from mrsimulator.utils import flatten_dict

from .base import ContribSchema

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


SITE_MAPPABLE = {
    "isotropic_chemical_shift": "isotropic",
    "shielding_symmetric.zeta": "zeta",
    "shielding_symmetric.eta": "eta",
    "shielding_symmetric.alpha": "alpha",
    "shielding_symmetric.beta": "beta",
    "shielding_symmetric.gamma": "gamma",
    "quadrupolar.Cq": "Cq",
    "quadrupolar.eta": "eta",
    "quadrupolar.alpha": "alpha",
    "quadrupolar.beta": "beta",
    "quadrupolar.gamma": "gamma",
}

NEW_MAPPABLE = {
    "shielding_symmetric": "ChemicalShift",
    "quadrupolar": "Quadrupolar",
    "isotropic_chemical_shift": "ChemicalShift",
}


def parse_sites(site):
    """Parse the site object and generate site contribution for mpcontribs.

    Args:
        Site site: Site object.
    """
    if site.shielding_symmetric is not None:
        site.shielding_symmetric.zeta *= -1  # convert to shift anisotropy
    site_ = flatten_dict(site.json())

    dict_ = {"ChemicalShift": {}, "Quadrupolar": {}}
    for k, v in site_.items():
        root = k.split(".")[0]
        new_key = SITE_MAPPABLE[k] if k in SITE_MAPPABLE else k
        if root in NEW_MAPPABLE:
            dict_[NEW_MAPPABLE[root]][new_key] = site_[k]
        else:
            dict_[new_key] = v

    _ = [dict_.pop(i) for i in ["ChemicalShift", "Quadrupolar"] if dict_[i] == {}]
    return dict_


def parse_method(method):
    """Parse the method object and generate method contribution for mpcontribs.

    Args:
        Method method: Method object.
    """
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


def mpcontribs_export(
    simulator: Simulator,
    signal_processors: list,
    project: str,
    identifier: str,
    exp_dict: dict = {},
    **kwargs,
):
    """Generate mpcontribs contribution entries for every site in the Simulator object.

    Arguments
    ---------
    Simulator simulator
        Object from where the site contributions are extracted.
    list signal_processors
        A list of SignalProcessor objects.
    str project
        mpcontribs project name (required).
    str identifier
        mpcontribs identifier (required).
    exp_dict
        Additional metadata to use in contribution.
    **kwargs
        Optional keyword arguments from mpcontribs ContributionsSchema.

    Example
    -------
    >>> contribution_data = mpcontribs_export(
    ...     simulator,
    ...     processors,
    ...     project='my_project',
    ...     identifier='mp-test'
    ... ) # doctest:+SKIP
    """

    attach = prepare_attachments(simulator, signal_processors, identifier)

    contribs = [
        ContribSchema(
            project=project,
            identifier=identifier,
            data={
                "site": {**parse_sites(site)},
                "method": {**parse_method(simulator.methods[0])},
            },
        ).json()
        for sys in simulator.spin_systems
        for site in sys.sites
    ]
    _ = [item.update({"attachments": attach}) for item in contribs]

    for item in contribs:
        item["data"]["method"] = {**item["data"]["method"], **exp_dict}

    return contribs


def save_obj(filename, data):
    with zopen(filename, "w") as f:
        json_str = json.dumps(data) + "\n"  # 2. string (i.e. JSON)
        json_bytes = json_str.encode("utf-8")
        f.write(json_bytes)


def load_obj(filename):
    with zopen(filename) as f:
        json_bytes = f.read()
        json_str = json_bytes.decode("utf-8")  # 2. string (i.e. JSON)
        data = json.loads(json_str)
        return data


def prepare_attachments(simulator, signal_processors, identifier):
    """Prepares the dict data as attachments for contribution."""
    filename = f"{identifier}.mrsim"
    mrsim = Attachment.from_data(filename, dict(simulator, signal_processors))
    data = [
        Attachment.from_data(f"{identifier}-exp{i}.csdf", mth.experiment.dict())
        for i, mth in enumerate(simulator.methods)
        if mth.experiment is not None
    ]
    return [mrsim, *data]
