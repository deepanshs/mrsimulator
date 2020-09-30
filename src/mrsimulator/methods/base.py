# -*- coding: utf-8 -*-
from copy import deepcopy

import numpy as np

from .utils import generate_method_from_template
from .utils import METHODS_DATA

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

# BlochDecaySpectrum
BlochDecaySpectrum = generate_method_from_template(METHODS_DATA["Bloch_decay"])

BlochDecayCentralTransitionSpectrum = generate_method_from_template(
    METHODS_DATA["Bloch_decay_central_transition"]
)


def Method2D(spectral_dimensions=[{}, {}], **kwargs):
    r"""A generic 2D correlation method.

    Attributes
    ----------

    channels: list (optional).
        The value is a list of isotope symbols over which the given method applies.
        An isotope symbol is given as a string with the atomic number followed by its
        atomic symbol, for example, '1H', '13C', and '33S'. The default is an empty
        list.
        The number of isotopes in a `channel` depends on the method. For example, a
        `BlochDecaySpectrum` method is a single channel method, in which case, the
        value of this attribute is a list with a single isotope symbol, ['13C'].

    rotor_angle: float (optional)
        A globally defined value for the angle between the sample rotation axis and the
        applied external magnetic field, :math:`\theta`, in units of rad. The default
        value is ``0.9553166``, i.e. the magic angle.

    magetic_flux_density: float (optional)
        A globally defined value for the macroscopic magnetic flux density, :math:`H_0`,
        of the applied external magnetic field in units of T. The default is ``9.4``.

    spectral_dimensions: list of :ref:`spectral_dim_api` or dict objects (optional).
        The number of spectral dimensions depends on the given method. For example, a
        `BlochDecaySpectrum` method is a one-dimensional method and thus requires a
        single spectral dimension. The default is a single default
        :ref:`spectral_dim_api` object.

    affine_matrix: np.ndarray or list (optional)
        An affine transformation square matrix,
        :math:`\mathbf{A} \in \mathbb{R}^{n \times n}`, where `n` is the number of
        spectral dimensions. The affine operation follows

        .. math::
            \mathbf{V}\prime = \mathbf{A} \cdot \mathbf{V},

        where :math:`\mathbf{V}\in\mathbb{R}^n` and :math:`\mathbf{V}\in\mathbb{R}^n`
        are the initial and transformed frequency coordinates.

    name: str (optional).
        The value is the name or id of the method. The default value is None.

    label: str (optional).
        The value is a label for the method. The default value is None.

    description: str (optional).
        The value is a description of the method. The default value is None.

    experiment: CSDM or ndarray (optional).
        An object holding the experimental measurement for the given method, if
        available. The default value is None.

    Return
    ------
        A :class:`~mrsimulator.Method` instance.

    Note
    ----

    The `rotor_frequency` parameter is fixed for this method and produces an
    infinite spinning speed spectrum.

    If the parameters `rotor_angle` and `magetic_flux_density` are defined outside
    of the `spectral_dimensions` list, the value of these parameters is considered
    global. In a multi-event method, you may also assign parameter values to
    individual events.
    """

    for dim in spectral_dimensions:
        if "events" in dim.keys():
            for evt in dim["events"]:
                if "transition_query" in evt.keys():
                    t_query = evt["transition_query"]
                    if "P" in t_query.keys():
                        if isinstance(t_query["P"], list):
                            t_query["P"] = {"channel-1": [[i] for i in t_query["P"]]}
                    if "D" in t_query.keys():
                        if isinstance(t_query["D"], list):
                            t_query["D"] = {"channel-1": [[i] for i in t_query["D"]]}

    Method2d_ = generate_method_from_template(deepcopy(METHODS_DATA["Method2D"]))
    return Method2d_(spectral_dimensions, **kwargs)


def message(attr, name):
    return f"`{attr}` attribute cannot be modified for {name} method."


def check_for_transition_query(name, spectral_dimensions=[{}, {}]):
    check = [
        "transition_query" in event.keys()
        for item in spectral_dimensions
        if "events" in item.keys()
        for event in item["events"]
    ]

    if np.any(check):
        raise AttributeError(message("transition_query", name))


def check_for_events(name, spectral_dimensions=[{}, {}]):
    check = ["events" in item.keys() for item in spectral_dimensions]

    if np.any(check):
        raise AttributeError(message("events", name))


# Args:
# channels: A list of isotope symbols over which the method will be applied.

# spectral_dimensions: A list of python dict. Each dict is contains keywords that
#     describe the coordinates along a spectral dimension. The keywords along with
#     its definition are:

#     - count:
#         An optional integer with the number of points, :math:`N`, along the
#         dimension. The default value is 1024.
#     - spectral_width:
#         An `optional` float with the spectral width, :math:`\Delta x`, along the
#         dimension in units of Hz. The default is 25 kHz.
#     - reference_offset:
#         An `optional` float with the reference offset, :math:`x_0` along the
#         dimension in units of Hz. The default value is 0 Hz.
#     - origin_offset:
#         An `optional` float with the origin offset (Larmor frequency) along the
#         dimension in units of Hz. The default value is None.
#     - events:
#         An `optional` list of Event objects. Each event object consists of
#         `magetic_flux_density`, `rotor_angle`, and `transition_query`
#         parameters, and are described below.

# rotor_angle: An `optional` float containing the angle between the sample
#     rotation axis and the applied external magnetic field, :math:`\theta`, in
#     units of rad. The default value is ``0.9553166``, i.e. the magic angle.
# magetic_flux_density: An `optional` float containing the macroscopic magnetic
#     flux density, :math:`H_0`, of the applied external magnetic field in units
#     of T. The default value is ``9.4``.
# affine_matrix: An `optional` affine matrix for affine transformation. The shape
#     of the affine matrix must be :math:`n \times n`, where `n` is the number of
#     spectral dimensions. The affine operation follows
