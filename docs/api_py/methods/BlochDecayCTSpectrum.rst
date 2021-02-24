
Bloch Decay Central Transition Spectrum method
----------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: BlochDecayCTSpectrum

**Example**

.. doctest::

    >>> from mrsimulator.methods import BlochDecayCTSpectrum
    >>> Bloch_CT_method = BlochDecayCTSpectrum(
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

Bloch decay central transition selective method is a special case of
:class:`~mrsimulator.methods.Method1D`, given as

.. doctest::

    >>> from mrsimulator.methods import Method1D
    >>> BlochdecayCT = BlochDecayCTSpectrum(
    ...     channels=['1H'],
    ...     rotor_frequency=5000, # in Hz
    ...     rotor_angle=0.95531, # in rad
    ...     magnetic_flux_density=9.4, # in T
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 50000, # in Hz
    ...             "reference_offset": -8000, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [-1], "D": [0]}
    ...             }]
    ...         }
    ...     ]
    ... )
