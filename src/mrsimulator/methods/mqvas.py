# -*- coding: utf-8 -*-
from copy import deepcopy

from .utils import generate_method_from_template
from .utils import METHODS_DATA

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class MQ_VAS_:
    r"""A generic multiple-quantum variable-angle spinning method for simulating average
    frequencies. The resulting spectrum is sheared such that the correlating dimensions
    are the isotropic dimension and the VAS dimension, respectively, where the isotropic
    dimension is given as

    .. math::
        \nu_\text{MQ-iso} = 1/(1+k) \nu_\text{MQ}(-m, m) + k/(1+k) \nu_\text{MAS}

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
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is ``0.9553166``, i.e. the magic angle.

    note:
        The `rotor_frequency` parameter is fixed for this method. The method produces
        an infinite spinning speed spectrum.

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
        template = deepcopy(METHODS_DATA["MQMAS_sheared"])
        template["name"] = cls.__name__
        method = generate_method_from_template(template)(spectral_dimensions, **kwargs)

        spin = method.channels[0].spin

        # select the coherence for the first event
        P = int(2 * mq)
        nQ = P
        P = -P if mq == spin else P

        method.spectral_dimensions[0].events[0].transition_query.P = {
            "channel-1": [[P]]
        }

        k = cls.k_MQ_MAS[nQ][spin]

        # Update the fractions for the events in the t1 spectral dimension.
        method.spectral_dimensions[0].events[0].fraction = 1 / (1 + k)
        method.spectral_dimensions[0].events[1].fraction = k / (1 + k)

        method.description = f"Simulate a {nQ}Q variable-angle spinning spectrum."
        return method


class ThreeQ_VAS(MQ_VAS_):
    r"""Simulate a sheared triple-quantum variable angle spinning spectrum correlating
    the symmetric triple-quantum transition,
    :math:`\left|\frac{3}{2}\right\rangle \rightarrow \left|-\frac{3}{2}\right\rangle`,
    to the central transition.

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
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is ``0.9553166``, i.e. the magic angle.

    note:
        The `rotor_frequency` parameter is fixed for this method. The method produces an
        infinite spinning speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> method = ThreeQ_VAS(
        ...     channels=["87Rb"],
        ...     magnetic_flux_density=7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 4e3,  # in Hz
        ...             "reference_offset": -5e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...         {
        ...             "count": 512,
        ...             "spectral_width": 1e4,  # in Hz
        ...             "reference_offset": -4e3,  # in Hz
        ...             "label": "MAS dimension",
        ...         },
        ...     ],
        ... )
    """

    def __new__(self, spectral_dimensions=[{}, {}], **kwargs):
        return super().__new__(
            self, mq=1.5, spectral_dimensions=spectral_dimensions, **kwargs
        )


class FiveQ_VAS(MQ_VAS_):
    def __new__(self, spectral_dimensions=[{}, {}], **kwargs):
        return super().__new__(
            self, mq=2.5, spectral_dimensions=spectral_dimensions, **kwargs
        )
