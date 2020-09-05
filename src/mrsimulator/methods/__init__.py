# -*- coding: utf-8 -*-
from copy import deepcopy

from .mqvas import FiveQ_VAS  # noqa:F401
from .mqvas import ThreeQ_VAS  # noqa:F401
from .stvas import ST1_VAS  # noqa:F401
from .utils import generate_method_from_template
from .utils import METHODS_DATA

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


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
