"""
Function list:
    - ThreeQ_VAS
    - FiveQ_VAS
    - SevenQ_VAS
"""
from mrsimulator.spin_system.isotope import Isotope

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

    note:
        The `rotor_frequency` parameter is fixed for this method. The method produces
        an infinite spinning speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.
    """

    @classmethod
    def update(cls, **kwargs):
        name = cls.__name__

        mq = MQ_p_symmetry[name]["mq"]
        iso = kwargs["channels"][0]
        spin = Isotope(symbol=iso).spin if isinstance(iso, str) else iso.spin

        # select the coherence for the first event
        P = int(2 * mq)
        nQ = P
        P = -P if mq == spin else P

        # set transition symmetry elements for spectral dimension 0
        events_0 = [{"transition_queries": [{"ch1": {"P": [P], "D": [0]}}]}]
        # set transition symmetry elements for spectral dimension 1
        events_1 = [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}]

        # method affine matrix shear factor
        k = shear_factor_MQ_MAS[nQ][spin]

        description = f"Simulate a {nQ}Q variable-angle spinning spectrum."
        return {
            "name": name,
            "description": description,
            "spectral_dimensions": [{"events": events_0}, {"events": events_1}],
            "affine_matrix": [1 / (1 + k), k / (1 + k), 0, 1],
        }


class ThreeQ_VAS(MQ_VAS):
    """Simulate a sheared and scaled 3Q 2D variable-angle spinning spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-1.5⟩⟨1.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    class Config:
        extra = "forbid"


class FiveQ_VAS(MQ_VAS):
    """Simulate a sheared and scaled 5Q variable-angle spinning spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-2.5⟩⟨2.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    class Config:
        extra = "forbid"


class SevenQ_VAS(MQ_VAS):
    """Simulate a sheared and scaled 7Q variable-angle spinning spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-3.5⟩⟨3.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    class Config:
        extra = "forbid"
