# -*- coding: utf-8 -*-
"""
Function list:
    - ST1_VAS
    - ST2_VAS
"""
from .base import BaseNamedMethod2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


shear_factor_ST_MAS = {
    3: {1.5: 24 / 27, 2.5: 21 / 72, 3.5: 84 / 135, 4.5: 165 / 216},
    5: {2.5: 132 / 72, 3.5: 69 / 135, 4.5: 12 / 216},
    7: {3.5: 324 / 135, 4.5: 243 / 216},
    9: {4.5: 600 / 216},
}

ST_p_symmetry = {"ST1_VAS": {"st": 1.5}, "ST2_VAS": {"st": 2.5}}


class ST_VAS(BaseNamedMethod2D):
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

    def __new__(cls, **kwargs):
        return cls.update(super().__new__(cls, **kwargs))

    @staticmethod
    def update(method):
        st = ST_p_symmetry[method.name]["st"]
        spin = method.channels[0].spin

        # select the coherence for the first event
        d = st ** 2 - (st - 1) ** 2
        D = [[d], [-d]]

        # setting transition symmetry elements
        sd = method.spectral_dimensions
        sd[0].events[0].transition_query.P = {"channel-1": [[-1]]}
        sd[0].events[0].transition_query.D = {"channel-1": D}

        sd[1].events[0].transition_query.P = {"channel-1": [[-1]]}
        sd[1].events[0].transition_query.D = {"channel-1": [[0]]}

        # method affine matrix
        if method.affine_matrix is None:
            k = shear_factor_ST_MAS[int(2 * st)][spin]
            method.affine_matrix = [1 / (1 + k), k / (1 + k), 0, 1]

        # method description
        method.description = (
            f"Simulate a {st} -> {st-1} and {-st+1} -> {-st} satellite-transition "
            "variable-angle spinning spectrum."
        )
        return method


class ST1_VAS(ST_VAS):
    r"""Simulate a sheared and scaled inner satellite and central transition correlation
    spectrum.

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
        >>> method = ST1_VAS(
        ...     channels=["87Rb"],
        ...     magnetic_flux_density=9.4,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 128,
        ...             "spectral_width": 1e3,  # in Hz
        ...             "reference_offset": -5e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 1e4,  # in Hz
        ...             "reference_offset": -3e3,  # in Hz
        ...             "label": "MAS dimension",
        ...         },
        ...     ],
        ... )
        >>> sys = SpinSystem(sites=[Site(isotope='87Rb')])
        >>> pprint(method.get_transition_pathways(sys))
        [|-1.5⟩⟨-0.5| ⟶ |-0.5⟩⟨0.5|, |0.5⟩⟨1.5| ⟶ |-0.5⟩⟨0.5|]
    """


class ST2_VAS(ST_VAS):
    r"""Simulate a sheared and scaled second to inner satellite and central transition
    correlation spectrum.

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
        >>> method = ST2_VAS(
        ...     channels=["17O"],
        ...     magnetic_flux_density=9.4,  # in T
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
        >>> sys = SpinSystem(sites=[Site(isotope='17O')])
        >>> pprint(method.get_transition_pathways(sys))
        [|-2.5⟩⟨-1.5| ⟶ |-0.5⟩⟨0.5|, |1.5⟩⟨2.5| ⟶ |-0.5⟩⟨0.5|]
    """
