
Bloch Decay Spectrum method
---------------------------

.. currentmodule:: mrsimulator.methods

.. autoclass:: BlochDecaySpectrum
    :show-inheritance:
    :members:
    :inherited-members: BaseModel

**Example**

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> from mrsimulator.method import SpectralDimension
    >>> Bloch_method = BlochDecaySpectrum(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         SpectralDimension(count=1024, spectral_width=50000, reference_offset=-8000)
    ...     ],
    ... )


Bloch decay method is a special case of :py:class:`~mrsimulator.methods.Method1D`, given as

.. doctest::

    >>> from mrsimulator.methods import Method1D
    >>> Blochdecay = Method1D(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 50000,  # in Hz
    ...             "reference_offset": -8000,  # in Hz
    ...             "events": [{"transition_query": [{"P": [-1]}]}],
    ...         }
    ...     ],
    ... )
