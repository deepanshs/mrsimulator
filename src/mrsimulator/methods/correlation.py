# -*- coding: utf-8 -*-
"""
Function list:
    - Cosy
    - Inadequate
"""
from .base import BaseNamedMethod2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class Cosy(BaseNamedMethod2D):
    r"""Simulate an infinite spinning COrrelation SpectroscopY spectrum.

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
            units of rad. The default value is 0.9553166, i.e. the magic angle.
    note:
        The `rotor_frequency` parameter is fixed for this method.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> cosy = Cosy(
        ...     channels=["13C"],
        ...     magnetic_flux_density=7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 4.4e3,  # in Hz
        ...             "reference_offset": 2.2e3,  # in Hz
        ...         },
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 4.4e3,  # in Hz
        ...             "reference_offset": 2.2e3,  # in Hz
        ...         },
        ...     ],
        ... )

    For two coupled spin-1/2 system, Cosy method selects 16 transition pathways.

        >>> sys = SpinSystem(
        ...     sites=[
        ...         Site(isotope='13C', isotropic_chemical_shift=10.0),
        ...         Site(isotope='13C', isotropic_chemical_shift=50.0,)
        ...     ],
        ...     couplings=[Coupling(site_index=[0, 1], isotropic_j=200)]
        ... )
        >>> pprint(cosy.get_transition_pathways(sys))
        [|-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|,
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|,
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|,
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|,
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|,
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|,
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|,
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|]
    """

    def __new__(cls, **kwargs):
        return cls.update(super().__new__(cls, **kwargs))

    @staticmethod
    def update(method):
        # setting transition symmetry elements
        sd = method.spectral_dimensions
        sd[0].events[0].transition_query.P = {"channel-1": [[-1]]}
        sd[1].events[0].transition_query.P = {"channel-1": [[-1]]}

        # method description
        method.description = (
            "Simulate an infinite spinning COrrelation SpectroscopY spectrum."
        )
        return method


class Inadequate(BaseNamedMethod2D):
    r"""Simulate an infinite spinning Incredible Natural Abundance DoublE QUAntum
    Transfer Experiment spectrum.

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
            units of rad. The default value is 0.9553166, i.e. the magic angle.
    note:
        The `rotor_frequency` parameter is fixed for this method.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> inadequate = Inadequate(
        ...     channels=["13C"],
        ...     magnetic_flux_density=7,  # in T
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 4.4e3,  # in Hz
        ...             "reference_offset": 2.2e3,  # in Hz
        ...         },
        ...         {
        ...             "count": 256,
        ...             "spectral_width": 4.4e3,  # in Hz
        ...             "reference_offset": 2.2e3,  # in Hz
        ...         },
        ...     ],
        ... )

    For two coupled spin-1/2 system, Cosy method selects 16 transition pathways.

        >>> sys = SpinSystem(
        ...     sites=[
        ...         Site(isotope='13C', isotropic_chemical_shift=10.0),
        ...         Site(isotope='13C', isotropic_chemical_shift=50.0,)
        ...     ],
        ...     couplings=[Coupling(site_index=[0, 1], isotropic_j=200)]
        ... )
        >>> pprint(inadequate.get_transition_pathways(sys))
        [|-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|,
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|]
    """

    def __new__(cls, **kwargs):
        return cls.update(super().__new__(cls, **kwargs))

    @staticmethod
    def update(method):
        # setting transition symmetry elements
        sd = method.spectral_dimensions
        sd[0].events[0].transition_query.P = {"channel-1": [[-1, -1]]}
        sd[1].events[0].transition_query.P = {"channel-1": [[-1]]}

        # method description
        method.description = (
            "Simulate an infinite spinning Incredible Natural Abundance DoublE QUAntum "
            "Transfer Experiment spectrum."
        )
        return method
