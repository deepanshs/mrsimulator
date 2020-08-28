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
    """Generate method object from json template."""
    # constructor
    def constructor(self, spectral_dimensions=[{}], parse=False, **kwargs):
        prep = prepare_method_structure(template, **kwargs)
        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        if common != set():
            info = "`, `".join(list(common))
            e = (
                f"The attribute(s) `{info}` cannot be modified for {prep['name']} "
                "method."
            )
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


class MQ_MAS_:
    r"""Simulate a triple-quantum magic-angle spinning spectrum.
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
        magetic_flux_density: An `optional` float containing the macroscopic magnetic
            flux density, :math:`H_0`, of the applied external magnetic field in units
            of T. The default value is ``9.4``.

    note:
        The `rotor_frequency` and `rotor_angle` parameters are fixed for this method.
        The method produces an infinite spinning speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.
    """

    k_MQ_MAS = {
        3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
        5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
        7: {3.5: 483 / 135, 4.5: 84 / 216},
        9: {4.5: 1116 / 216},
    }

    def __new__(cls, mq=1.5, spectral_dimensions=[{}, {}], **kwargs):
        if "rotor_angle" in kwargs:
            e = "`rotor_angle` is fixed to the magic-angle and cannot be modified."
            raise AttributeError(e)

        template = deepcopy(METHODS_DATA["MQMAS_sheared"])
        template["name"] = cls.__name__
        method = generate_method_from_template(template)(spectral_dimensions, **kwargs)

        spin = method.channels[0].spin

        # select the coherence for the first event
        p = int(2 * mq)
        P = -p if mq == spin else p

        method.spectral_dimensions[0].events[0].transition_query.P["channel-1"] = [[P]]
        method.spectral_dimensions[0].events[0].transition_query.D["channel-1"] = [[0]]

        k = cls.k_MQ_MAS[p][spin]

        # Update the fractions for the events in the t1 spectral dimension.
        method.spectral_dimensions[0].events[0].fraction = 1 / (1 + k)
        method.spectral_dimensions[0].events[1].fraction = k / (1 + k)

        method.description = f"Simulate a {p}Q magic-angle spinning spectrum."
        return method


class ThreeQ_MAS(MQ_MAS_):
    def __new__(self, spectral_dimensions=[{}, {}], **kwargs):
        return super().__new__(
            self, mq=1.5, spectral_dimensions=spectral_dimensions, **kwargs
        )


class FiveQ_MAS(MQ_MAS_):
    def __new__(self, spectral_dimensions=[{}, {}], **kwargs):
        return super().__new__(
            self, mq=2.5, spectral_dimensions=spectral_dimensions, **kwargs
        )
