"""
Function list:
    - ST1_VAS
    - ST2_VAS
"""
from mrsimulator.spin_system.isotope import Isotope

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
    """Simulate a satellite-transition magic-angle spinning spectrum.

    note:
        The `rotor_frequency` and `rotor_angle` parameters are fixed for this method.
        The method produces an infinite spinning speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.
    """

    class Config:
        extra = "forbid"

    @classmethod
    def update(cls, **kwargs):
        name = cls.__name__
        st = ST_p_symmetry[name]["st"]
        description = (
            f"Simulate a {st} -> {st-1} and {-st+1} -> {-st} satellite-transition "
            "variable-angle spinning spectrum."
        )
        spin = Isotope(symbol=kwargs["channels"][0]).spin

        # select the coherence for the first event
        d = st**2 - (st - 1) ** 2

        # setting transition symmetry elements for spectral dimension 0
        events_0 = [
            {
                "transition_queries": [
                    {"ch1": {"P": [-1], "D": [d]}},
                    {"ch1": {"P": [-1], "D": [-d]}},
                ]
            }
        ]
        # setting transition symmetry elements for spectral dimension 1
        events_1 = [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}]

        # method affine matrix shear factor
        k = shear_factor_ST_MAS[int(2 * st)][spin]
        sign = 1  # if st == spin else -1

        return {
            "name": name,
            "description": description,
            "spectral_dimensions": [{"events": events_0}, {"events": events_1}],
            "affine_matrix": [1 / (1 + k), sign * k / (1 + k), 0, 1],
        }


class ST1_VAS(ST_VAS):
    """Simulate a sheared and scaled inner satellite and central transition correlation
    spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-1.5⟩⟨-0.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j),
         |0.5⟩⟨1.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    class Config:
        extra = "forbid"


class ST2_VAS(ST_VAS):
    """Simulate a sheared and scaled second to inner satellite and central transition
    correlation spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-2.5⟩⟨-1.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j),
         |1.5⟩⟨2.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    class Config:
        extra = "forbid"
