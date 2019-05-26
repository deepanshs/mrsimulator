

.. _spectrum:

--------
Spectrum
--------

In version 0.1, the `Spectrum` object is a python
`dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object with a single key, `direct_dimension`, whose value is a
:ref:`direct_dimension` object.


*Example of Spectrum object.*


    >>> spectrum_object = {
    ...     "direct_dimension": {
    ...         "nucleus": "13C",
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "5 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "100 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }
