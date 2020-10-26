.. _generic_2d_method:

Generic two-dimensional correlation method
------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: Method2D

**Example**

.. doctest::

    >>> from mrsimulator.methods import Method2D
    >>> method = Method2D(
    ...     channels=["87Rb"],
    ...     magnetic_flux_density=7,  # in T. Global value for `magnetic_flux_density`.
    ...     rotor_angle=0.95531, # in rads. Global value for the `rotor_angle`.
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 256,
    ...             "spectral_width": 4e3,  # in Hz
    ...             "reference_offset": -5e3,  # in Hz
    ...             "event": [
    ...                 { # Global value for the `magnetic_flux_density` and `rotor_angle`
    ...                   # is used during this event.
    ...                     "transition_query": {"P": [-3], "D": [0]}
    ...                 }
    ...             ]
    ...         },
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 1e4,  # in Hz
    ...             "reference_offset": -4e3,  # in Hz
    ...             "event": [
    ...                 { # Global value for the `magnetic_flux_density` is used during this
    ...                   # event. User defined local value for `rotor_angle` is used here.
    ...                     "rotor_angle": 1.2238, # in rads
    ...                     "transition_query": {"P": [-1], "D": [0]}
    ...                 }
    ...             ]
    ...         },
    ...     ],
    ...     affine_matrix=[[1, -1], [0, 1]],
    ... )
