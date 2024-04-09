from typing import List

from .base import BaseNamedMethod2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class SSB2D(BaseNamedMethod2D):
    r"""Simulating a sheared 2D finite to infinite speed MAS correlation spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-0.5⟩⟨0.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]
    """

    name: str = "SSB2D"
    description: str = "Simulate a 2D sideband separation method."
    affine_matrix: List = [1, -1, 0, 1]

    class Config:
        extra = "forbid"

    def __init__(self, **kwargs):
        error = f"`rotor_frequency` cannot be zero for {__class__.__name__} method."
        if "rotor_frequency" not in kwargs:
            raise ValueError(error)
        if kwargs["rotor_frequency"] == 0:
            raise ValueError(error)

        super().__init__(**kwargs)

    @classmethod
    def update(cls, **kwargs):
        # events for spectral dimension at index 0
        events_0 = [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}]
        # events for spectral dimension at index 1
        events_1 = [
            {
                "rotor_frequency": 1.0e12,
                "transition_queries": [{"ch1": {"P": [-1], "D": [0]}}],
            },
        ]
        return {"spectral_dimensions": [{"events": events_0}, {"events": events_1}]}
