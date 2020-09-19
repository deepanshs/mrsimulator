# -*- coding: utf-8 -*-
from copy import deepcopy
from os import path

from monty.serialization import loadfn
from mrsimulator.method import Method
from mrsimulator.method.event import Event
from mrsimulator.method.spectral_dimension import SpectralDimension

MODULE_DIR = path.dirname(path.abspath(__file__))
METHODS_DATA = loadfn(path.join(MODULE_DIR, "methods_data.json"))

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

# create class dynamically
__doc_args__ = r"""

Args:
    spectral_dimensions: A list of python dict. Each dict is contains keywords that
        describe the coordinates along a spectral dimension. The keywords along with
        its definition are:

        count:
            An optional integer with the number of points, :math:`N`, along the
            dimension. The default value is 1024.
        spectral_width:
            An `optional` float with the spectral width, :math:`\Delta x`, along the
            dimension in units of Hz. The default is 25 kHz.
        reference_offset:
            An `optional` float with the reference offset, :math:`x_0` along the
            dimension in units of Hz. The default value is 0 Hz.
        origin_offset:
            An `optional` float with the origin offset (Larmor frequency) along the
            dimension in units of Hz. The default value is None.

    channels: A list of isotope symbols over which the method will be applied.
    rotor_frequency: An `optional` float containing the sample spinning frequency
        :math:`\nu_r`, in units of Hz. The default value is ``0``.
    rotor_angle: An `optional` float containing the angle between the sample rotation
        axis and the applied external magnetic field, :math:`\theta`, in units of rad.
        The default value is ``0.9553166``, i.e. the magic angle.
    magetic_flux_density: An `optional` float containing the macroscopic magnetic
        flux density, :math:`H_0`, of the applied external magnetic field in units of
        T. The default value is ``9.4``.

Return:
    A :class:`~mrsimulator.Method` instance.
"""


def prepare_method_structure(template, **kwargs):
    keys = kwargs.keys()
    n_channels = template["number_of_channels"]
    name = template["name"] if "name" not in keys else kwargs["name"]
    desc = (
        template["description"] if "description" not in keys else kwargs["description"]
    )
    label = None if "label" not in keys else kwargs["label"]
    experiment = None if "experiment" not in keys else kwargs["experiment"]
    affine_matrix = None if "affine_matrix" not in keys else kwargs["affine_matrix"]

    prep = {
        "name": name,
        "description": desc,
        "label": label,
        "experiment": experiment,
        "affine_matrix": affine_matrix,
    }

    if "channels" in kwargs:
        prep["channels"] = kwargs["channels"]
        given_n_channels = len(prep["channels"])
        if given_n_channels != n_channels:
            raise ValueError(
                f"The method requires exactly {n_channels} channel(s), "
                f"{given_n_channels} provided."
            )
        kwargs.pop("channels")
    else:
        prep["channels"] = ["1H" for _ in range(n_channels)]
    return prep


def generate_method_from_template(template):
    """Generate method object from json template."""
    # constructor
    def constructor(self, spectral_dimensions=[{}], **kwargs):
        parse = False
        if "parse" in kwargs:
            parse = kwargs["parse"]
        prep = prepare_method_structure(template, **kwargs)
        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        if common != set():
            info = "`, `".join(list(common))
            e = f"`{info}` attribute cannot be modified for {prep['name']} class."
            raise AttributeError(e)

        dim = []
        n_sp = len(spectral_dimensions)
        n_tem = len(template["spectral_dimensions"])

        if n_tem < n_sp:
            raise ValueError(
                f"The method allows {n_tem} spectral dimension(s), {n_sp} given."
            )

        for i, s in enumerate(template["spectral_dimensions"]):
            events = []

            _fill_missing_events_in_template(spectral_dimensions[i], s)

            for j, e in enumerate(s["events"]):
                ew = set(e)

                kw = kwargs.copy()
                if "events" in spectral_dimensions[i]:
                    for key, val in spectral_dimensions[i]["events"][j].items():
                        kw[key] = val

                intersection = ew.intersection(kw)
                [e.pop(item) for item in intersection]

                # prioritize the keyword arguments over the global arguments.
                common = set(kw).intersection(set(global_events))
                ge = deepcopy(global_events)
                [ge.pop(item) for item in common]

                params = {**e, **kw, **ge}
                params = params if parse else Event(**params)
                events.append(params)

            if "events" in spectral_dimensions[i]:
                spectral_dimensions[i].pop("events")

            params = {**spectral_dimensions[i], "events": events}
            params = params if parse else SpectralDimension(**params)
            dim.append(params)

        method = {**prep, "spectral_dimensions": dim}
        method = Method.parse_dict_with_units(method) if parse else Method(**method)
        return method

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of the method object from as a python dictionary.

        Args:
            py_dict: Dict object.
        """
        dict_copy = deepcopy(py_dict)
        spectral_dimensions = dict_copy["spectral_dimensions"]
        dict_copy.pop("spectral_dimensions")
        return cls.__new__(
            0, spectral_dimensions=spectral_dimensions, parse=True, **dict_copy
        )

    description = template["description"]
    method = type(
        template["name"],
        (object,),
        {
            "__new__": constructor,
            "__str__": description,
            "__doc__": description + __doc_args__,
            "parse_dict_with_units": parse_dict_with_units,
        },
    )
    return method


def _fill_missing_events_in_template(spectral_dimensions, s_template):
    """Fill the missing events in the template relative to the spectral dimensions."""
    if "events" not in spectral_dimensions:
        return

    s_tem_len = len(s_template["events"])
    sp_evt_len = len(spectral_dimensions["events"])
    if s_tem_len < sp_evt_len:
        diff = sp_evt_len - s_tem_len
        [s_template["events"].append(Event().dict()) for _ in range(diff)]
