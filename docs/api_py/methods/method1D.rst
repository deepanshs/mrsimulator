.. _generic_1d_method:

Generic one-dimensional method
------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: Method1D

**Example**

Simulating triple-quantum 1D spectrum.

.. doctest::

    >>> from mrsimulator.methods import Method1D
    >>> method1 =  Method1D(
    ...     channels=["87Rb"],
    ...     magnetic_flux_density=7,  # in T
    ...     rotor_angle=54.735*np.pi/180,
    ...     rotor_frequency=1e9,
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 1e4,  # in Hz
    ...             "reference_offset": -4e3,  # in Hz
    ...             "label": "quad only",
    ...             "events":[
    ...                 {"transition_query": {"P": [-3], "D": [0]}},
    ...             ]
    ...         }
    ...     ]
    ... )
