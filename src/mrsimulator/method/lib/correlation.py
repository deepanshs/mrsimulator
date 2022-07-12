"""
Function list:
    - Cosy
    - Inadequate
"""
import numpy as np

from .base import BaseNamedMethod2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class Cosy(BaseNamedMethod2D):
    """Simulate an infinite spinning COrrelation SpectroscopY spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|, weight=(0.25+0j),
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(-0.25-0j),
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(0.25+0j),
         |-0.5, -0.5⟩⟨-0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|, weight=(0.25+0j),
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|, weight=(-0.25+0j),
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(0.25+0j),
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(0.25+0j),
         |-0.5, -0.5⟩⟨0.5, -0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|, weight=(0.25+0j),
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|, weight=(0.25+0j),
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(0.25+0j),
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(0.25+0j),
         |-0.5, 0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|, weight=(-0.25+0j),
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|, weight=(0.25+0j),
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(0.25+0j),
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(-0.25+0j),
         |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|, weight=(0.25+0j)]
    """

    name: str = "Cosy"
    description: str = (
        "Simulate an infinite spinning COrrelation SpectroscopY spectrum."
    )

    @classmethod
    def update(self, **kwargs):
        event_0 = [
            {"transition_queries": [{"ch1": {"P": [-1]}}]},
            {"query": {"ch1": {"angle": np.pi / 2}}},
        ]
        event_1 = [{"transition_queries": [{"ch1": {"P": [-1]}}]}]

        return {
            "spectral_dimensions": [{"events": event_0}, {"events": event_1}],
        }


class Inadequate(BaseNamedMethod2D):
    """Simulate an infinite spinning Incredible Natural Abundance DoublE QUAntum
    Transfer Experiment spectrum.

    Note:
        The attribute `rotor_frequency` cannot be modified for this method and is set to
        simulate an infinite speed spectrum.

    Return:
        A :py:class:`~mrsimulator.Method` instance.

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
        [|-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨-0.5, 0.5|, weight=(1+0j),
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(1+0j),
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(1+0j),
         |-0.5, -0.5⟩⟨0.5, 0.5| ⟶ |0.5, -0.5⟩⟨0.5, 0.5|, weight=(1+0j)]
    """

    name: str = "Inadequate"
    description: str = (
        "Simulate an infinite spinning Incredible Natural Abundance DoublE QUAntum "
        "Transfer Experiment spectrum."
    )

    @classmethod
    def update(cls, **kwargs):
        event_0 = [{"transition_queries": [{"ch1": {"P": [-1, -1]}}]}]
        event_1 = [{"transition_queries": [{"ch1": {"P": [-1]}}]}]

        return {
            "spectral_dimensions": [{"events": event_0}, {"events": event_1}],
        }
