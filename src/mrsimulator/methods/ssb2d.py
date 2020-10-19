# -*- coding: utf-8 -*-
from . import base as bs

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def SSB2D(**kwargs):
    r"""A specialized method for simulating 2D finite speed to infinite speed MAS
    correlation spectum. For spin I=1/2, the infinite speed MAS is the isotropic
    dimension. The resulting spectrum is sheared.

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

        rotor_frequency: An `optional` float containing the sample spinning frequency
            :math:`\nu_r`, in units of Hz. The default value is ``0``.
        magetic_flux_density: An `optional` float containing the macroscopic magnetic
            flux density, :math:`H_0`, of the applied external magnetic field in units
            of T. The default value is 9.4.
        rotor_angle: An `optional` float containing the angle between the sample
            rotation axis and the applied external magnetic field, :math:`\theta`, in
            units of rad. The default value is 0.9553166, i.e. the magic angle.

    Return:
        A :class:`~mrsimulator.Method` instance.

    Example:
        >>> method = SSB2D(
        ...     channels=["13C"],
        ...     magnetic_flux_density=7,  # in T
        ...     rotor_frequency=1500, # in Hz
        ...     spectral_dimensions=[
        ...         {
        ...             "count": 16,
        ...             "spectral_width": 16*1500,  # in Hz (= count * rotor_frequency)
        ...             "reference_offset": -5e3,  # in Hz
        ...             "label": "Sideband dimension",
        ...         },
        ...         {
        ...             "count": 512,
        ...             "spectral_width": 1e4,  # in Hz
        ...             "reference_offset": -4e3,  # in Hz
        ...             "label": "Isotropic dimension",
        ...         },
        ...     ],
        ... )
        >>> sys = SpinSystem(sites=[Site(isotope='13C')])
        >>> method.get_transition_pathways(sys)
        [TransitionPathway(|-0.5⟩⟨0.5|, |-0.5⟩⟨0.5|)]
    """
    name = "SSB2D"

    spectral_dimensions = bs.check_for_spectral_dimensions(kwargs, 2)
    bs.check_for_events(name, spectral_dimensions)

    # check for rotor frequency
    error = f"`rotor_frequency` cannot be zero for {name} method."
    if "rotor_frequency" not in kwargs:
        raise ValueError(error)
    if kwargs["rotor_frequency"] == 0:
        raise ValueError(error)

    spin_freq = kwargs["rotor_frequency"]
    kwargs.pop("rotor_frequency")

    method = bs.Method2D(spectral_dimensions, name=name, **kwargs)
    method.spectral_dimensions[0].events[0].rotor_frequency = spin_freq
    method.spectral_dimensions[1].events[0].rotor_frequency = 1e9

    # check that the spectral width for the first dimension is equal to
    # spin frequency * count
    # count = method.spectral_dimensions[0].count
    # if method.spectral_dimensions[0].spectral_width != spin_freq * count:
    #     raise ValueError(
    #         "The spectral width along the dimension at index 0 must be equal to the "
    #         "product of the `spin_frequency` and `count` of the dimension."
    #     )

    method.description = "Simulate a 2D sideband separation method."

    # Add the affine matrix
    if method.affine_matrix is None:
        method.affine_matrix = [1, -1, 0, 1]

    return method
