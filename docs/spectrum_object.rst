

.. _spectrum:

---------------
Spectrum object
---------------

In version 0.1, the `Spectrum` is a python dict object with a single key,
``direct_dimension``, whose value is a `DirectDimension` object.


**DirectDimension object**
---

  A `DirectDimension` object is a python dict object with the following
  key-value pairs.

  .. list-table::
    :widths: 25 75
    :header-rows: 1

    * - key
      - value description
    * - ``nucleus``
      - The value is a string with the isotope symbol of the nuclei
        offset of the spectrum, e.g. '1.4587 kHz'.

    * - ``magnetic_flux_density``
      - The value is a string with a physical quantity describing the strength
        of the magnetic filed of the NMR spectrometer, e.g. '9.4 T'.
    * - ``rotor_frequency``
      - The value is a string with a physical quantity describing the sample
        spinning frequency, e.g. '10 kHz'.
    * - ``rotor_angle``
      - The value is a string with a physical quantity describing the angle
        of the sample rotation axis with respect to the external magnetic field,
        e.g. '54.735 deg'.
    * - ``number_of_points``
      - The value is an interger with the number of points used for sampling the
        NMR spectrum, e.g. 8192.
    * - ``spectral_width``
      - The value is a string with a physical quantity containing the frequency
        spectral width of the spectrum, e.g. '100 kHz'.
    * - ``reference_offset``
      - The value is a string with a physical quantity containing the reference
        offset of the spectrum, e.g. '1.4587 kHz'.

    * - ``rotor_phase``
      - The value is a a string with a physical quantity describing the angle
        of the sample rotation axis with respect to the external magnetic field.



**Example of Spectrum object**

    >>> {
    ...     "direct_dimension": {
    ...         "nucleus": "13C"
    ...         "magnetic_flux_density": "9.4 T",
    ...         "rotor_frequency": "5 kHz",
    ...         "rotor_angle": "54.735 deg",
    ...         "number_of_points": 8192,
    ...         "spectral_width": "100 kHz",
    ...         "reference_offset": "0 Hz",
    ...     }
    ... }
