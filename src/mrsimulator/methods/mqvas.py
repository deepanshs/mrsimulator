# -*- coding: utf-8 -*-
"""
Function list:
    - ThreeQ_VAS
    - FiveQ_VAS
    - SevenQ_VAS
"""
from .base import BaseNamedMethod2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


shear_factor_MQ_MAS = {
    3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
    5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
    7: {3.5: 483 / 135, 4.5: 84 / 216},
    9: {4.5: 1116 / 216},
}

MQ_p_symmetry = {
    "ThreeQ_VAS": {"mq": 1.5},
    "FiveQ_VAS": {"mq": 2.5},
    "SevenQ_VAS": {"mq": 3.5},
}


class MQ_VAS(BaseNamedMethod2D):
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
            of T. The default value is 9.4.
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is 0.9553166, i.e. the magic angle.

    note:
        The `rotor_frequency` parameter is fixed for this method. The method produces
        an infinite spinning speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.
    """

    def __new__(cls, **kwargs):
        return cls.update(super().__new__(cls, **kwargs))

    @staticmethod
    def update(method):
        mq = MQ_p_symmetry[method.name]["mq"]
        spin = method.channels[0].spin

        # select the coherence for the first event
        P = int(2 * mq)
        nQ = P
        P = -P if mq == spin else P

        # setting transition symmetry elements
        sd = method.spectral_dimensions
        sd[0].events[0].transition_query.P = {"channel-1": [[P]]}
        sd[0].events[0].transition_query.D = {"channel-1": [[0]]}

        sd[1].events[0].transition_query.P = {"channel-1": [[-1]]}
        sd[1].events[0].transition_query.D = {"channel-1": [[0]]}

        # method affine matrix
        if method.affine_matrix is None:
            k = shear_factor_MQ_MAS[nQ][spin]
            method.affine_matrix = [1 / (1 + k), k / (1 + k), 0, 1]

        # method description
        method.description = f"Simulate a {nQ}Q variable-angle spinning spectrum."
        return method


class ThreeQ_VAS(MQ_VAS):
    r"""Simulate a sheared and scaled 3Q 2D variable-angle spinning spectrum.

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
            of T. The default value is 9.4.
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is 0.9553166, i.e. the magic angle.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> method = ThreeQ_VAS(
        ...     channels=["87Rb"],
        ...     magnetic_flux_density=7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 128,
        ...             "spectral_width": 3e3,  # in Hz
        ...             "reference_offset": -2e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...         {
        ...             "count": 512,
        ...             "spectral_width": 1e4,  # in Hz
        ...             "reference_offset": -5e3,  # in Hz
        ...             "label": "MAS dimension",
        ...         },
        ...     ],
        ... )
        >>> sys = SpinSystem(sites=[Site(isotope='87Rb')])
        >>> method.get_transition_pathways(sys)
        [|-1.5⟩⟨1.5| ⟶ |-0.5⟩⟨0.5|]
    """


class FiveQ_VAS(MQ_VAS):
    r"""Simulate a sheared and scaled 5Q variable-angle spinning spectrum.

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
            of T. The default value is 9.4.
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is 0.9553166, i.e. the magic angle.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> method = FiveQ_VAS(
        ...     channels=["17O"],
        ...     magnetic_flux_density=11.7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 512,
        ...             "spectral_width": 5e3,  # in Hz
        ...             "reference_offset": -3e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...         {
        ...             "count": 512,
        ...             "spectral_width": 2e4,  # in Hz
        ...             "reference_offset": -2e3,  # in Hz
        ...             "label": "MAS dimension",
        ...         },
        ...     ],
        ... )
        >>> sys = SpinSystem(sites=[Site(isotope='17O')])
        >>> method.get_transition_pathways(sys)
        [|-2.5⟩⟨2.5| ⟶ |-0.5⟩⟨0.5|]
    """


class SevenQ_VAS(MQ_VAS):
    r"""Simulate a sheared and scaled 7Q variable-angle spinning spectrum.

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
            of T. The default value is 9.4.
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is 0.9553166, i.e. the magic angle.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> method = SevenQ_VAS(
        ...     channels=["51V"],
        ...     magnetic_flux_density=4.7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 1e3,  # in Hz
        ...             "reference_offset": 1e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...         {
        ...             "count": 1024,
        ...             "spectral_width": 1e4,  # in Hz
        ...             "reference_offset": -2e3,  # in Hz
        ...             "label": "MAS dimension",
        ...         },
        ...     ],
        ... )
        >>> sys = SpinSystem(sites=[Site(isotope='51V')])
        >>> method.get_transition_pathways(sys)
        [|-3.5⟩⟨3.5| ⟶ |-0.5⟩⟨0.5|]
    """
