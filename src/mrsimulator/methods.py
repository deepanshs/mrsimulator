# -*- coding: utf-8 -*-
from copy import deepcopy
from os import path

from monty.serialization import loadfn

from .method import Event
from .method import Method
from .method import SpectralDimension

# program to create class dynamically

MODULE_DIR = path.dirname(path.abspath(__file__))
METHODS_DATA = loadfn(path.join(MODULE_DIR, "methods_data.json"))

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

    channels: A list of isotope symbols over which the method will be applied. The
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


def generate_method_from_template(template, __doc__):
    # constructor
    def constructor(self, spectral_dimensions=[{}], parse=False, **kwargs):
        n_channels = template["number_of_channels"]
        prep = {"name": template["name"], "description": template["description"]}
        if "channels" in kwargs:
            prep["channels"] = kwargs["channels"]
            given_n_channels = len(prep["channels"])
            if given_n_channels != n_channels:
                raise ValueError(
                    f"The method requires exactly {n_channels} channels, "
                    f"{given_n_channels} provided."
                )
            kwargs.pop("channels")
        else:
            prep["channels"] = ["1H" for _ in range(n_channels)]

        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        if common != set():
            raise AttributeError(
                f"The attribute(s) {common} is fixed for {prep['name']} method."
            )

        dim = []
        for i, s in enumerate(template["spectral_dimensions"]):
            events = []
            for e in s["events"]:
                ew = set(e)
                intersection = ew.intersection(kw)
                kw = kwargs.copy()
                for item in intersection:
                    kw.pop(item)
                if parse:
                    events.append({**e, **kw, **global_events})
                else:
                    events.append(Event(**e, **kw, **global_events))
            if parse:
                dim.append({**spectral_dimensions[i], "events": events})
            else:
                dim.append(SpectralDimension(**spectral_dimensions[i], events=events))

        if parse:
            method = {**prep, "spectral_dimensions": dim}
            return Method.parse_dict_with_units(method)
        else:
            return Method(**prep, spectral_dimensions=dim)

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

    method = type(
        template["name"],
        (object,),
        {
            "__new__": constructor,
            "__str__": template["description"],
            "__doc__": __doc__ + __doc_args__,
            "parse_dict_with_units": parse_dict_with_units,
        },
    )
    return method


BlochDecaySpectrum = generate_method_from_template(
    METHODS_DATA["Bloch_decay"], "Method for simulating Bloch decay spectrum."
)

BlochDecayCentralTransitionSpectrum = generate_method_from_template(
    METHODS_DATA["Bloch_decay_central_transition"],
    "Method for simulating central transition selective Bloch decay spectrum.",
)
