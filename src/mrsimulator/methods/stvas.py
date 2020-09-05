# -*- coding: utf-8 -*-
from copy import deepcopy

from .utils import generate_method_from_template
from .utils import METHODS_DATA

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class ST_VAS_:
    r"""Simulate a satellite-transition magic-angle spinning spectrum.
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

    k_ST_MAS = {
        3: {1.5: 24 / 27, 2.5: 21 / 72, 3.5: 84 / 135, 4.5: 165 / 216},
        5: {2.5: 132 / 72, 3.5: 69 / 135, 4.5: 12 / 216},
        7: {3.5: 324 / 135, 4.5: 243 / 216},
        9: {4.5: 600 / 216},
    }

    def __new__(cls, st=1.5, spectral_dimensions=[{}, {}], **kwargs):
        # if "rotor_angle" in kwargs:
        #     e = "`rotor_angle` is fixed to the magic-angle and cannot be modified."
        #     raise AttributeError(e)

        template = deepcopy(METHODS_DATA["MQMAS_sheared"])
        template["name"] = cls.__name__
        method = generate_method_from_template(template)(spectral_dimensions, **kwargs)

        spin = method.channels[0].spin

        # select the coherence for the first event
        d = st ** 2 - (st - 1) ** 2
        D = [[d], [-d]]

        method.spectral_dimensions[0].events[0].transition_query.D = {"channel-1": D}

        k = cls.k_ST_MAS[int(2 * st)][spin]

        # Update the fractions for the events in the t1 spectral dimension.
        method.spectral_dimensions[0].events[0].fraction = 1 / (1 + k)
        method.spectral_dimensions[0].events[1].fraction = k / (1 + k)

        method.description = (
            f"Simulate a {st} -> {st-1} and {-st+1} -> {-st} satellite-transition "
            "magic-angle spinning spectrum."
        )
        return method


class ST1_VAS(ST_VAS_):
    """Simulate a inner satellite to central transition ST-MAS correlation spectrum.

    The inner satellite are |3/2> -> |1/2> and |-1/2> -> |-3/2>.
    """

    def __new__(self, spectral_dimensions=[{}, {}], **kwargs):
        return super().__new__(
            self, st=1.5, spectral_dimensions=spectral_dimensions, **kwargs
        )
