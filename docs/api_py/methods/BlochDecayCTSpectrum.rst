
Bloch Decay Central Transition Spectrum method
----------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: BlochDecayCentralTransitionSpectrum

**Example**

.. doctest::

    >>> from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
    >>> Bloch_CT_method = BlochDecayCentralTransitionSpectrum(
    ...     channels=['1H'],
    ...     rotor_frequency=5000, # in Hz
    ...     rotor_angle=0.95531, # in rad
    ...     magnetic_flux_density=9.4, # in T
    ...     spectral_dimensions=[dict(
    ...         count=1024,
    ...         spectral_width=50000, # in Hz
    ...         reference_offset=-8000, # in Hz
    ...     )]
    ... )
