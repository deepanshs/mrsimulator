
Bloch Decay Central Transition Spectrum method
----------------------------------------------

.. currentmodule:: mrsimulator.method.lib

.. autoclass:: BlochDecayCTSpectrum
    :show-inheritance:
    :members:
    :inherited-members: BaseModel

**Example**

.. doctest::

    >>> from mrsimulator.method.lib import BlochDecayCTSpectrum
    >>> from mrsimulator.method import SpectralDimension
    >>> bloch_decay_ct = BlochDecayCTSpectrum(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         SpectralDimension(count=1024, spectral_width=50000, reference_offset=-8000)
    ...     ],
    ... )

Bloch decay central transition selective method is part of the methods library,
but can be constructed with a generic method as follows

.. doctest::

    >>> from mrsimulator.method import Method
    >>> bloch_decay_ct = Method(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 50000,  # in Hz
    ...             "reference_offset": -8000,  # in Hz
    ...             "events": [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}],
    ...         }
    ...     ],
    ... )
