
Bloch Decay Central Transition Spectrum method
----------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autoclass:: BlochDecayCTSpectrum
    :show-inheritance:
    :members:
    :inherited-members: BaseModel

**Example**

.. doctest::

    >>> from mrsimulator.methods import BlochDecayCTSpectrum
    >>> from mrsimulator.method import SpectralDimension
    >>> Bloch_CT_method = BlochDecayCTSpectrum(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         SpectralDimension(count=1024, spectral_width=50000, reference_offset=-8000)
    ...     ],
    ... )

Bloch decay central transition selective method is a special case of
:py:class:`~mrsimulator.methods.Method1D`, given as

.. doctest::

    >>> from mrsimulator.methods import Method1D
    >>> BlochdecayCT = Method1D(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 50000,  # in Hz
    ...             "reference_offset": -8000,  # in Hz
    ...             "events": [{"transition_query": [{"P": [-1], "D": [0]}]}],
    ...         }
    ...     ],
    ... )
