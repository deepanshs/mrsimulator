# -*- coding: utf-8 -*-
from copy import deepcopy
from os import path

from monty.serialization import loadfn
from mrsimulator.method import Method
from mrsimulator.method.event import Event
from mrsimulator.method.spectral_dimension import SpectralDimension

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

# create class dynamically

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
    n_channels = template["number_of_channels"]
    prep = {"name": template["name"], "description": template["description"]}
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
    # constructor
    def constructor(self, spectral_dimensions=[{}], parse=False, **kwargs):
        prep = prepare_method_structure(template, **kwargs)
        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        e = f"The attribute(s) {common} connot be modified for {prep['name']} method."
        if common != set():
            raise AttributeError(e)

        dim = []
        for i, s in enumerate(template["spectral_dimensions"]):
            events = []
            for j, e in enumerate(s["events"]):
                ew = set(e)

                kw = kwargs.copy()
                if "events" in spectral_dimensions[i]:
                    for key, val in spectral_dimensions[i]["events"][j].items():
                        kw[key] = val
                    spectral_dimensions[i].pop("events")

                intersection = ew.intersection(kw)
                [e.pop(item) for item in intersection]

                params = {**e, **kw, **global_events}
                params = params if parse else Event(**params)
                events.append(params)

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


# BlochDecaySpectrum
BlochDecaySpectrum = generate_method_from_template(METHODS_DATA["Bloch_decay"])

BlochDecayCentralTransitionSpectrum = generate_method_from_template(
    METHODS_DATA["Bloch_decay_central_transition"]
)

Custom2D = generate_method_from_template(METHODS_DATA["Custom2D"])


def MQVAS(spectral_dimensions=[{}, {}], **kwargs):
    r"""Simulate a multi-quantum variable angle spinning spectrum.

    Args:
        channels: A list of isotope symbols over which the method will be applied.

        spectral_dimensions: A list of python dict. Each dict is contains keywords that
            describe the coordinates along a spectral dimension. The keywords along with
            its definition are:

            - count:
                An optional integer with the number of points, :math:`N`, along the
                dimension. The default value is 1024.
            - spectral_width:
                An `optional` float with the spectral width, :math:`\Delta x`, along the
                dimension in units of Hz. The default is 25 kHz.
            - reference_offset:
                An `optional` float with the reference offset, :math:`x_0` along the
                dimension in units of Hz. The default value is 0 Hz.
            - origin_offset:
                An `optional` float with the origin offset (Larmor frequency) along the
                dimension in units of Hz. The default value is None.
            - events:
                An `optional` list of Event objects. Each event object consists of
                `magetic_flux_density`, `rotor_angle`, and `transition_query`
                parameters, and are described below.

        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is ``0.9553166``, i.e. the magic angle.
        magetic_flux_density: An `optional` float containing the macroscopic magnetic
            flux density, :math:`H_0`, of the applied external magnetic field in units
            of T. The default value is ``9.4``.

    note:
        The `rotor_frequency` parameter is fixed for this method and produces an
        infinite spinning speed spectrum.

        If the parameters `rotor_angle` and `magetic_flux_density` are defined outside
        of the `spectral_dimensions` list, the value of these parameters is considered
        global. In a multi-event method, such as the two-dimensional methods, you can
        also assign parameter values to individual events.

    Return:
        A :class:`~mrsimulator.Method` instance.

    """

    for dim in spectral_dimensions:
        if "events" in dim.keys():
            for evt in dim["events"]:
                if "transition_query" in evt.keys():
                    t_query = evt["transition_query"]
                    if "P" in t_query.keys():
                        t_query["P"] = {"channel-1": [[i] for i in t_query["P"]]}
                    if "D" in t_query.keys():
                        t_query["D"] = {"channel-1": [[i] for i in t_query["D"]]}

    MQVAS_ = generate_method_from_template(deepcopy(METHODS_DATA["MQVAS"]))
    return MQVAS_(spectral_dimensions, **kwargs)
