
Bloch Decay Spectrum method
---------------------------

.. currentmodule:: mrsimulator.method.lib

.. autoclass:: BlochDecaySpectrum
    :show-inheritance:
    :members:
    :inherited-members: BaseModel

**Example**

.. doctest::

    >>> from mrsimulator.method.lib import BlochDecaySpectrum
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


Bloch decay method is is part of the methods library,
but can be constructed with a generic method as follows

.. doctest::

    >>> from mrsimulator.method import Method
    >>> Blochdecay = Method(
    ...     channels=["1H"],
    ...     rotor_frequency=5000,  # in Hz
    ...     rotor_angle=54.735 * 3.14159 / 180,  # in rad
    ...     magnetic_flux_density=9.4,  # in T
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 50000,  # in Hz
    ...             "reference_offset": -8000,  # in Hz
    ...             "events": [{"transition_query": [{"ch1": {"P": [-1]}}]}],
    ...         }
    ...     ],
    ... )
